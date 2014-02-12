//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \file DTK_AsynchronousMap_def.hpp
 * \author Stuart R. Slattery
 * \brief Shared domain map definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ASYNCHRONOUSMAP_DEF_HPP
#define DTK_ASYNCHRONOUSMAP_DEF_HPP

#include <algorithm>
#include <limits>

#include "DTK_FieldTools.hpp"
#include "DTK_DBC.hpp"
#include "DTK_MeshTools.hpp"
#include "DTK_AKDTree.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Ptr.hpp>

#include <Tpetra_Distributor.hpp>
#include <Tpetra_MultiVector.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param comm The communicator over which the map is generated.
 *
 * \param dimension The dimension of the map. This should be consistent with
 * all source and target objects (i.e. only 3 dimensional coordinates will be
 * accepted with a 3 dimensional map). We need this here so we have a global
 * baseline for all objects that may or may not exist on all processes.
 *
 * \param store_missed_points Set to true if it is desired to keep track of
 * the local target points missed during map generation. The default value is
 * false. 
 */
template<class Mesh, class CoordinateField, int DIM>
AsynchronousMap<Mesh,CoordinateField,DIM>::AsynchronousMap( 
    const RCP_Comm& comm, bool store_missed_points )
    : d_comm( comm )
    , d_store_missed_points( store_missed_points )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Mesh, class CoordinateField, int DIM>
AsynchronousMap<Mesh,CoordinateField,DIM>::~AsynchronousMap()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Generate the shared domain map.
 *
 * \param source_mesh_manager Source mesh in the shared domain problem. A null
 * RCP is a valid argument. This will be the case when a mesh manager is only
 * constructed on a subset of the processes that the shared domain map is
 * constructed over. Note that the source mesh must exist only on processes
 * that reside within the AsynchronousMap communicator.
 *
 * \param target_coord_manager Target coordinates in the shared domain
 * problem. A null RCP is a valid argument. This will be the case when a field
 * manager is only constructed on a subset of the processes that the shared
 * domain map is constructed over. Note that the target coordinates must exist
 * only on processes that reside within the AsynchronousMap communicator.
 *
 * \param tolerance Absolute tolerance for point searching. Will be used when
 * checking the reference cell ( and is therefore absolute ).
 */
template<class Mesh, class CoordinateField, int DIM>
void AsynchronousMap<Mesh,CoordinateField,DIM>::setup( 
    const RCP_MeshManager& source_mesh_manager, 
    const RCP_CoordFieldManager& target_coord_manager,
    const int max_buffer_size,
    const int buffer_check_frequency,
    double tolerance )
{
    // Create existence values for the managers.
    bool source_exists = Teuchos::nonnull( source_mesh_manager );
    bool target_exists = Teuchos::nonnull( target_coord_manager );

    // Create local to global process indexer for the source.
    RCP_Comm source_comm;
    if ( source_exists )
    {
	source_comm = source_mesh_manager->comm();
    }
    d_source_indexer = CommIndexer( d_comm, source_comm );

    // Build global indexing for the target points.
    buildTargetMap( target_coord_manager, target_exists );
    Teuchos::ArrayView<const GlobalOrdinal> target_ordinals =
	d_target_map->getNodeElementList();

    // Create the communication plan.
    Teuchos::Array<int> source_neighbor_ranks;
    Teuchos::Array<BoundingBox> source_neighbor_boxes;
    Teuchos::Array<int> target_neighbor_ranks;
    createCommunicationPlan( source_mesh_manager,
			     source_exists,
			     target_coord_manager,
			     target_exists,
			     source_neighbor_ranks,
			     source_neighbor_boxes,
			     target_neighbor_ranks );

    // Create the asynchronous tree and perform point location.
    Teuchos::Array<GlobalOrdinal> target_point_gids;
    {
	AKDTree<Mesh,CoordinateField,DIM> akd_tree( 
	    d_comm, max_buffer_size, buffer_check_frequency,
	    source_neighbor_ranks, source_neighbor_boxes,
	    target_neighbor_ranks,
	    source_mesh_manager, target_coord_manager,
	    target_ordinals );

	akd_tree.locate( d_source_elements,
			 target_point_gids,
			 d_target_coords,
			 d_target_map->getGlobalNumElements(),
			 tolerance );
    }

    // Build the source map from the target ordinals.
    Teuchos::ArrayView<const GlobalOrdinal> target_point_gids_view = 
	target_point_gids();
    d_source_map = Tpetra::createNonContigMap<int,GlobalOrdinal>( 
	target_point_gids_view, d_comm );
    DTK_ENSURE( Teuchos::nonnull(d_source_map) );

    // Build the source-to-target exporter.
    d_source_to_target_exporter = 
      Teuchos::rcp( new Tpetra::Export<int,GlobalOrdinal>(
			  d_source_map, d_target_map ) );
    DTK_ENSURE( Teuchos::nonnull(d_source_to_target_exporter) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the points missed in the map generation.
 *
 * \return If store_missed_points is true, return the local indices of the
 *  points provided by target_coord_manager that were not mapped. An exception
 *  will be thrown if store_missed_points is false. Returns a null view if all
 *  points have been mapped or the map has not yet been generated.
*/
template<class Mesh, class CoordinateField, int DIM>
Teuchos::ArrayView<const typename 
		   AsynchronousMap<Mesh,CoordinateField,DIM>::GlobalOrdinal> 
AsynchronousMap<Mesh,CoordinateField,DIM>::getMissedTargetPoints() const
{
    DTK_REQUIRE( d_store_missed_points );
    
    return d_missed_points();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the points missed in the map generation.
 *
 * \return If store_missed_points is true, return the local indices of the
 *  points provided by target_coord_manager that were not mapped. An exception
 *  will be thrown if store_missed_points is false. Returns a null view if all
 *  points have been mapped or the map has not yet been generated.
*/
template<class Mesh, class CoordinateField, int DIM>
Teuchos::ArrayView<typename 
		   AsynchronousMap<Mesh,CoordinateField,DIM>::GlobalOrdinal> 
AsynchronousMap<Mesh,CoordinateField,DIM>::getMissedTargetPoints()
{
    DTK_REQUIRE( d_store_missed_points );
    
    return d_missed_points();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply the shared domain map for a valid source field evaluator and
 * target data space to the target points that were mapped.
 *
 * \param source_evaluator Function evaluator used to apply the mapping. This
 * FieldEvaluator must be valid for the source mesh used to generate the map.
 *
 * \param target_space_manager Target space into which the function
 * evaluations will be written. Enough space must be allocated to hold
 * evaluations at all points in all dimensions of the field.
 */
template<class Mesh, class CoordinateField, int DIM>
template<class SourceField, class TargetField>
void AsynchronousMap<Mesh,CoordinateField,DIM>::apply( 
    const Teuchos::RCP< FieldEvaluator<GlobalOrdinal,SourceField> >& source_evaluator,
    Teuchos::RCP< FieldManager<TargetField> >& target_space_manager )
{
    typedef FieldTraits<SourceField> SFT;
    typedef FieldTraits<TargetField> TFT;

    // Set existence values for the source and target.
    bool source_exists = Teuchos::nonnull( source_evaluator );
    bool target_exists = Teuchos::nonnull( target_space_manager );

    // Evaluate the source function at the target points and construct a view
    // of the function evaluations.
    int field_dim = 0;
    Teuchos::ArrayRCP<typename SFT::value_type> source_field_copy(0,0);
    if ( source_exists )
    {
	SourceField function_evaluations = 
	    source_evaluator->evaluate( 
		Teuchos::arcpFromArray( d_source_elements ),
		Teuchos::arcpFromArray( d_target_coords ) );

	field_dim = SFT::dim( function_evaluations );

	source_field_copy =    
	    FieldTools<SourceField>::copy( function_evaluations );
    }
    Teuchos::broadcast<int,int>( *d_comm, d_source_indexer.l2g(0),
				 Teuchos::Ptr<int>(&field_dim) );

    // Build a multivector for the function evaluations.
    GlobalOrdinal source_size = source_field_copy.size() / field_dim;
    Teuchos::RCP<Tpetra::MultiVector<typename SFT::value_type, int, GlobalOrdinal> > 
	source_vector = Tpetra::createMultiVectorFromView( 
	    d_source_map, source_field_copy, source_size, field_dim );

    // Construct a view of the target space. Fill the target space with zeros
    // so that points we didn't map get some data.
    Teuchos::ArrayRCP<typename TFT::value_type> target_field_view(0,0);
    GlobalOrdinal target_size = 0;
    if ( target_exists )
    {
	target_field_view = FieldTools<TargetField>::nonConstView( 
	    *target_space_manager->field() );

	target_size = target_field_view.size() / field_dim;

	FieldTools<TargetField>::putScalar( 
	    *target_space_manager->field(), 0.0 );

	DTK_CHECK( 
	    target_size == Teuchos::as<GlobalOrdinal>(
		d_target_map->getNodeNumElements()) );
    }

    // Build a multivector for the target space.
    Teuchos::RCP<Tpetra::MultiVector<typename TFT::value_type, int, GlobalOrdinal> > 
	target_vector =	Tpetra::createMultiVectorFromView( 
	    d_target_map, target_field_view, target_size, field_dim );

    // Move the data from the source decomposition to the target
    // decomposition.
    target_vector->doExport( *source_vector, *d_source_to_target_exporter, 
			     Tpetra::INSERT );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the target map from the input target coordinates. Each
 * evaluation point is assumed unique.
 */
template<class Mesh, class CoordinateField, int DIM>
void AsynchronousMap<Mesh,CoordinateField,DIM>::buildTargetMap(
    const RCP_CoordFieldManager& target_coord_manager,
    bool target_exists )
{
    // Get the local number of target points.
    GlobalOrdinal local_num_targets = 0;
    if ( target_exists )
    {
	local_num_targets = FieldTools<CoordinateField>::dimSize( 
	    *target_coord_manager->field() );
    }

    // Get the total number of target points in the entire problem.
    GlobalOrdinal global_num_targets = 0;
    Teuchos::reduceAll<int,GlobalOrdinal>( *d_comm,
					   Teuchos::REDUCE_SUM,
					   1,
					   &local_num_targets,
					   &global_num_targets );

    // Build the data import map from the point global ordinals.
    d_target_map = Tpetra::createContigMap<int,GlobalOrdinal>(
	global_num_targets, local_num_targets, d_comm );
    DTK_ENSURE( Teuchos::nonnull(d_target_map) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the target map from the input target coordinates. Each
 * evaluation point is assumed unique.
 */
template<class Mesh, class CoordinateField, int DIM>
void AsynchronousMap<Mesh,CoordinateField,DIM>::createCommunicationPlan(
    const RCP_MeshManager& source_mesh_manager,
    const bool source_exists,
    const RCP_CoordFieldManager& target_coord_manager,
    const bool target_exists,
    Teuchos::Array<int>& source_neighbor_ranks,
    Teuchos::Array<BoundingBox>& source_neighbor_boxes,
    Teuchos::Array<int>& target_neighbor_ranks )
{
    // Build the local source bounding box.
    BoundingBox local_source_box;
    if ( source_exists )
    {
	local_source_box = source_mesh_manager->localBoundingBox();
    }

    // Gather the source boxes.
    Teuchos::Array<BoundingBox> source_boxes( d_comm->getSize() );
    Teuchos::gatherAll<int,BoundingBox>( *d_comm,
					 1,
					 &local_source_box,
					 source_boxes.size(),
					 source_boxes.getRawPtr() );

    // Build the local target bounding box and find the source procs that
    // neighbor this target proc. These are our neighbors that we will
    // send buffers to.
    if ( target_exists )
    {
	BoundingBox local_target_box
	    = FieldTools<CoordinateField>::coordLocalBoundingBox(
		*target_coord_manager->field() );

	for ( unsigned i = 0; i < source_boxes.size(); ++i )
	{
	    if ( BoundingBox::checkForIntersection(
		     local_target_box,source_boxes[i]) )
	    {
		source_neighbor_ranks.push_back(i);
		source_neighbor_boxes.push_back(source_boxes[i]);
	    }
	}
    }

    // Build the communication plan so the source procs know the
    // neighboring target procs they will receive buffers from.
    Tpetra::Distributor distributor( d_comm );
    distributor.createFromSends( source_neighbor_ranks );
    target_neighbor_ranks = distributor.getImagesFrom();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ASYNCHRONOUSMAP_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_AsynchronousMap_def.hpp
//---------------------------------------------------------------------------//
