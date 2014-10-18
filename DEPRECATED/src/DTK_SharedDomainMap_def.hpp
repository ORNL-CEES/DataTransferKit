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
 * \file DTK_SharedDomainMap_def.hpp
 * \author Stuart R. Slattery
 * \brief Shared domain map definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SHAREDDOMAINMAP_DEF_HPP
#define DTK_SHAREDDOMAINMAP_DEF_HPP

#include <algorithm>
#include <limits>

#include "DTK_FieldTools.hpp"
#include "DTK_DBC.hpp"
#include "DTK_Rendezvous.hpp"
#include "DTK_MeshTools.hpp"

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
template<class CoordinateField>
SharedDomainMap<CoordinateField>::SharedDomainMap( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm, 
    const int dimension, 
    bool store_missed_points )
    : d_comm( comm )
    , d_dimension( dimension )
    , d_store_missed_points( store_missed_points )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class CoordinateField>
SharedDomainMap<CoordinateField>::~SharedDomainMap()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Generate the shared domain map.
 *
 * \param source_mesh_manager Source mesh in the shared domain problem. A null
 * RCP is a valid argument. This will be the case when a mesh manager is only
 * constructed on a subset of the processes that the shared domain map is
 * constructed over. Note that the source mesh must exist only on processes
 * that reside within the SharedDomainMap communicator.
 *
 * \param target_coord_manager Target coordinates in the shared domain
 * problem. A null RCP is a valid argument. This will be the case when a field
 * manager is only constructed on a subset of the processes that the shared
 * domain map is constructed over. Note that the target coordinates must exist
 * only on processes that reside within the SharedDomainMap communicator.
 *
 * \param tolerance Reference tolerance for point searching. Will be used when
 * checking the reference cell for point inclusion.
 */
template<class CoordinateField>
void SharedDomainMap<CoordinateField>::setup( 
    const Teuchos::RCP<MeshManager>& source_mesh_manager, 
    const RCP_CoordFieldManager& target_coord_manager,
    double tolerance )
{
    // Create existence values for the managers.
    bool source_exists = Teuchos::nonnull( source_mesh_manager );
    bool target_exists = Teuchos::nonnull( target_coord_manager );

    // Create local to global process indexers for the managers.
    Teuchos::RCP<const Teuchos::Comm<int> > source_comm;
    if ( source_exists )
    {
	source_comm = source_mesh_manager->comm();
    }
    d_source_indexer = CommIndexer( d_comm, source_comm );
    Teuchos::RCP<const Teuchos::Comm<int> > target_comm;
    if ( target_exists )
    {
	target_comm = target_coord_manager->comm();
    }
    d_target_indexer = CommIndexer( d_comm, target_comm );

    // Create a communicator for asynchronous communication.
    Teuchos::RCP<const Teuchos::Comm<int> > async_comm = d_comm->duplicate();

    // Check the source and target dimensions for consistency and build the
    // global bounding boxes.
    BoundingBox source_box;
    if ( source_exists )
    {
	DTK_REQUIRE( source_mesh_manager->dim() == d_dimension );
	source_box = source_mesh_manager->globalBoundingBox();
    }
    BoundingBox local_target_box( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
    BoundingBox send_target_box;
    BoundingBox receive_target_box;
    if ( target_exists )
    {
	DTK_REQUIRE( CFT::dim( *target_coord_manager->field() ) 
		     == d_dimension );
	send_target_box = FieldTools<CoordinateField>::coordGlobalBoundingBox(
	    *target_coord_manager->field(), target_coord_manager->comm() );
	local_target_box = FieldTools<CoordinateField>::coordLocalBoundingBox(
	    *target_coord_manager->field() );
    }

    // Post a receive for the target box on source proc 0.
    Teuchos::RCP<Teuchos::CommRequest<int> > box_request;
    if ( async_comm->getRank() == d_source_indexer.l2g(0) )
    {
	box_request = Teuchos::ireceive<int,BoundingBox>(
	    *async_comm, Teuchos::RCP<BoundingBox>(&receive_target_box,false),
	    d_target_indexer.l2g(0) );
    }

    // Send the target box to source proc 0 from target proc 0.
    if ( async_comm->getRank() == d_target_indexer.l2g(0) )
    {
	Teuchos::isend<int,BoundingBox>( 
	    *async_comm, Teuchos::RCP<BoundingBox>(&send_target_box,false), 
	    d_source_indexer.l2g(0) );
    }

    // Intersect the boxes on source proc 0 to get the shared domain bounding
    // box.
    BoundingBox shared_domain_box;
    if ( async_comm->getRank() == d_source_indexer.l2g(0) )
    {
        Teuchos::Ptr<Teuchos::RCP<Teuchos::CommRequest<int> > > 
	    request_ptr(&box_request);
        Teuchos::wait( *async_comm, request_ptr );
	bool has_intersect = BoundingBox::intersectBoxes( 
	    source_box, receive_target_box, shared_domain_box );
	DTK_INSIST( has_intersect );
    }
    Teuchos::broadcast<int,BoundingBox>( 
	*d_comm, d_source_indexer.l2g(0),
	Teuchos::Ptr<BoundingBox>(&shared_domain_box) );

    // Build a rendezvous decomposition with the source mesh.
    Rendezvous rendezvous( d_comm, d_dimension, shared_domain_box );
    rendezvous.build( source_mesh_manager, d_source_indexer, local_target_box );

    // Get the local number of target points.
    MeshId local_num_targets = 0;
    if ( target_exists )
    {
	local_num_targets = FieldTools<CoordinateField>::dimSize( 
	    *target_coord_manager->field() );
    }

    // Get the total number of target points in the entire problem.
    MeshId global_num_targets = 0;
    Teuchos::reduceAll<int,MeshId>( *d_comm,
				    Teuchos::REDUCE_SUM,
				    1,
				    &local_num_targets,
				    &global_num_targets );

    // Build the data import map from the point global ordinals.
    d_target_map = Tpetra::createContigMap<int,MeshId>(
	global_num_targets, local_num_targets, d_comm );
    DTK_ENSURE( Teuchos::nonnull(d_target_map) );

    // Get the target ordinals from the map.
    Teuchos::ArrayView<const MeshId> target_ordinals =
	d_target_map->getNodeElementList();
    DTK_CHECK( local_num_targets == target_ordinals.size() );

    // If we're keeping track of missed points, we also need to build the
    // global-to-local ordinal map.
    if ( d_store_missed_points )
    {
	for ( MeshId n = 0; n < local_num_targets; ++n )
	{
	    d_target_g2l[ target_ordinals[n] ] = n;
	}
    }

    // Determine the rendezvous destination proc of each point in the
    // coordinate field. Also get the target points that are in the box in
    // which the rendezvous decomposition was generated. The rendezvous
    // algorithm will expand the box slightly based on mesh parameters.
    Teuchos::Array<double> coords_copy(0,0.0);
    Teuchos::Array<MeshId> targets_in_box;
    if ( target_exists )
    {
	coords_copy.resize( CFT::size(*target_coord_manager->field()) );
	std::copy( CFT::begin(*target_coord_manager->field()),
		   CFT::end(*target_coord_manager->field()),
		   coords_copy.begin() );

	getTargetPointsInBox( rendezvous.getBox(), 
			      *target_coord_manager->field(),
			      target_ordinals, targets_in_box );
    }
    Teuchos::Array<int> rendezvous_procs = 
	rendezvous.procsContainingPoints( Teuchos::arcpFromArray(coords_copy) );

    // Extract those target points that are not in the box. We don't want to
    // send these to the rendezvous decomposition.
    Teuchos::Array<MeshId> not_in_box;
    for ( unsigned i = 0; i < targets_in_box.size(); ++i )
    {
	if ( targets_in_box[i] == std::numeric_limits<MeshId>::max() )
	{
	    not_in_box.push_back( i );
	}
    }
    std::reverse( not_in_box.begin(), not_in_box.end() );

    typename Teuchos::Array<MeshId>::const_iterator not_in_box_iterator;
    unsigned remove_stride = 0;
    for ( not_in_box_iterator = not_in_box.begin();
	  not_in_box_iterator != not_in_box.end();
	  ++not_in_box_iterator )
    {
	remove_stride = rendezvous_procs.size();

	rendezvous_procs.remove( *not_in_box_iterator );

	for ( int d = d_dimension - 1; d >= 0; --d )
	{
	    coords_copy.remove( 
		d*remove_stride + (*not_in_box_iterator) );
	}
    }
    not_in_box.clear();
    DTK_CHECK( coords_copy.size() % d_dimension == 0 );
    DTK_CHECK( rendezvous_procs.size() ==
	       coords_copy.size() / d_dimension );

    typename Teuchos::Array<MeshId>::iterator targets_bound =
	std::remove( targets_in_box.begin(), targets_in_box.end(), 
		     std::numeric_limits<MeshId>::max() );
    MeshId targets_in_box_size = 
	std::distance( targets_in_box.begin(), targets_bound );

    targets_in_box.resize( targets_in_box_size );

    // Via an inverse communication operation, move the global point ordinals
    // that are in the rendezvous decomposition box to the rendezvous
    // decomposition.
    Teuchos::ArrayView<const MeshId> targets_in_box_view = 
	targets_in_box();
    Tpetra::Distributor target_to_rendezvous_distributor( d_comm );
    MeshId num_rendezvous_points = 
	target_to_rendezvous_distributor.createFromSends( rendezvous_procs() );
    Teuchos::Array<MeshId> rendezvous_points( num_rendezvous_points );
    target_to_rendezvous_distributor.doPostsAndWaits( 
	targets_in_box_view, 1, rendezvous_points() );

    // Move the target coordinates to the rendezvous decomposition.
    Teuchos::Array<double> 
	rendezvous_coords( d_dimension*num_rendezvous_points );
    MeshId num_points = rendezvous_procs.size();
    Teuchos::ArrayView<const double> coords_dim;
    Teuchos::ArrayView<double> points_dim;
    for ( int d = 0; d < d_dimension; ++d )
    {
	coords_dim = coords_copy( d*num_points, num_points );
	points_dim = rendezvous_coords( d*num_rendezvous_points, 
					num_rendezvous_points );
	target_to_rendezvous_distributor.doPostsAndWaits(
	    coords_dim, 1, points_dim );
    }
    coords_copy.clear();

    // Search the rendezvous decomposition with the target points to get the
    // source elements that contain them.
    Teuchos::Array<MeshId> rendezvous_elements;
    Teuchos::Array<int> rendezvous_element_src_procs;
    rendezvous.elementsContainingPoints( 
	Teuchos::arcpFromArray(rendezvous_coords),
	rendezvous_elements,
	rendezvous_element_src_procs,
	tolerance );

    // Get the points that were not in the mesh. If we're keeping track of
    // missed points, also make a list of those ordinals.
    Teuchos::Array<MeshId> not_in_mesh;
    Teuchos::Array<MeshId> missed_in_mesh_idx;
    Teuchos::Array<MeshId> missed_in_mesh_ordinal;
    for ( unsigned i = 0; i < rendezvous_elements.size(); ++i )
    {
	if ( rendezvous_elements[i] == 
	     std::numeric_limits<MeshId>::max() )
	{
	    not_in_mesh.push_back( i );

	    if ( d_store_missed_points )
	    {
		missed_in_mesh_idx.push_back( i );
		missed_in_mesh_ordinal.push_back( rendezvous_points[i] );
	    }
	}
    }

    // If we're keeping track of missed points, send their global ordinals
    // back to the target decomposition so that we can add them to the list.
    if ( d_store_missed_points )
    {
	// Extract the missed point target procs from the target-to-rendezvous
	// distributor.
	Teuchos::ArrayView<const int> from_images = 
	    target_to_rendezvous_distributor.getImagesFrom();
	Teuchos::ArrayView<const std::size_t> from_lengths = 
	    target_to_rendezvous_distributor.getLengthsFrom();
	Teuchos::Array<int> point_target_procs;
	for ( int i = 0; i < (int) from_images.size(); ++i )
	{
	    for ( std::size_t j = 0; j < from_lengths[i]; ++j )
	    {
		point_target_procs.push_back( from_images[i] );
	    }
	}
	DTK_CHECK( Teuchos::as<MeshId>(point_target_procs.size())
		   == num_rendezvous_points );

	// Build a list of target procs for the missed points.
	Teuchos::Array<int> missed_target_procs( missed_in_mesh_idx.size() );
	for ( int n = 0; n < (int) missed_in_mesh_idx.size(); ++n )
	{
	    missed_target_procs[n] = 
		point_target_procs[ missed_in_mesh_idx[n] ];
	}
	point_target_procs.clear();

	// Send the missed points back to the target decomposition through an
	// inverse communication operation and add them to the list.
	Teuchos::ArrayView<const MeshId> missed_in_mesh_ordinal_view = 
	    missed_in_mesh_ordinal();
	Tpetra::Distributor rendezvous_to_target_distributor( d_comm );
	MeshId num_missed_targets = 
	    rendezvous_to_target_distributor.createFromSends( 
		missed_target_procs() );
	MeshId offset = d_missed_points.size();
	d_missed_points.resize( offset + num_missed_targets );
	rendezvous_to_target_distributor.doPostsAndWaits( 
	    missed_in_mesh_ordinal_view, 1, 
	    d_missed_points.view( offset, num_missed_targets ) );

	// Convert the missed point global indices to local indices and add
	// them to the list.
	for ( MeshId n = offset; n < offset+num_missed_targets; ++n )
	{
	    d_missed_points[n] = 
		d_target_g2l.find( d_missed_points[n] )->second;
	}
    }
    missed_in_mesh_idx.clear();
    missed_in_mesh_ordinal.clear();

    // Extract the points we didn't find in any elements in the rendezvous
    // decomposition and their corresponding elements. We don't want to send
    // these to the source.
    std::reverse( not_in_mesh.begin(), not_in_mesh.end() );
    typename Teuchos::Array<MeshId>::const_iterator not_in_mesh_iterator;
    for ( not_in_mesh_iterator = not_in_mesh.begin();
	  not_in_mesh_iterator != not_in_mesh.end();
	  ++not_in_mesh_iterator )
    {
	remove_stride = rendezvous_points.size();
	rendezvous_points.remove( *not_in_mesh_iterator );
	for ( int d = d_dimension - 1; d >= 0; --d )
	{
	    rendezvous_coords.remove( 
		d*remove_stride + (*not_in_mesh_iterator) );
	}
    }
    not_in_mesh.clear();
    DTK_CHECK( rendezvous_coords.size() % d_dimension == 0 );
    DTK_CHECK( rendezvous_points.size() ==
	       rendezvous_coords.size() / d_dimension );

    typename Teuchos::Array<MeshId>::iterator rendezvous_elements_bound =
	std::remove( rendezvous_elements.begin(), 
		     rendezvous_elements.end(), 
		     std::numeric_limits<MeshId>::max() );
    MeshId rendezvous_elements_size = 
	std::distance( rendezvous_elements.begin(), rendezvous_elements_bound );

    typename Teuchos::Array<int>::iterator rendezvous_element_src_procs_bound =
	std::remove( rendezvous_element_src_procs.begin(), 
		     rendezvous_element_src_procs.end(), -1 );
    MeshId rendezvous_element_src_procs_size = 
	std::distance( rendezvous_element_src_procs.begin(), 
		       rendezvous_element_src_procs_bound );
    DTK_CHECK( rendezvous_elements_size == 
	       rendezvous_element_src_procs_size );

    rendezvous_elements.resize( rendezvous_elements_size );
    rendezvous_element_src_procs.resize( rendezvous_element_src_procs_size );

    // Setup rendezvous-to-source distributor.
    Tpetra::Distributor rendezvous_to_src_distributor( d_comm );
    MeshId num_source_elements = 
	rendezvous_to_src_distributor.createFromSends( 
	    rendezvous_element_src_procs() );

    // Send the rendezvous elements to the source decomposition via inverse
    // communication. 
    Teuchos::ArrayView<const MeshId> rendezvous_elements_view =
	rendezvous_elements();
    d_source_elements.resize( num_source_elements );
    rendezvous_to_src_distributor.doPostsAndWaits( 
	rendezvous_elements_view, 1, d_source_elements() );

    // Send the rendezvous point global ordinals to the source decomposition
    // via inverse communication.
    Teuchos::ArrayView<const MeshId> reduced_rendezvous_points_view =
	rendezvous_points();
    Teuchos::Array<MeshId> source_points( num_source_elements );
    rendezvous_to_src_distributor.doPostsAndWaits( 
	reduced_rendezvous_points_view, 1, source_points() );

    // Build the source map from the target ordinals.
    Teuchos::ArrayView<const MeshId> source_points_view = 
	source_points();
    d_source_map = Tpetra::createNonContigMap<int,MeshId>( 
	source_points_view, d_comm );
    DTK_ENSURE( Teuchos::nonnull(d_source_map) );

    // Send the rendezvous point coordinates to the source decomposition.
    d_target_coords.resize( num_source_elements*d_dimension );
    Teuchos::ArrayView<const double> const_points_dim;
    Teuchos::ArrayView<double> target_coord_dim;
    num_rendezvous_points = rendezvous_points.size();
    for ( int d = 0; d < d_dimension; ++d )
    {
	const_points_dim = rendezvous_coords( d*num_rendezvous_points, 
					      num_rendezvous_points );
	target_coord_dim = d_target_coords( d*num_source_elements,
					    num_source_elements );
	rendezvous_to_src_distributor.doPostsAndWaits(
	    const_points_dim, 1, target_coord_dim );
    }

    // Build the source-to-target exporter.
    d_source_to_target_exporter = 
	Teuchos::rcp( new Tpetra::Export<int,MeshId>(
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
template<class CoordinateField>
Teuchos::ArrayView<const MeshId> 
SharedDomainMap<CoordinateField>::getMissedTargetPoints() const
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
template<class CoordinateField>
Teuchos::ArrayView<MeshId> 
SharedDomainMap<CoordinateField>::getMissedTargetPoints()
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
template<class CoordinateField>
template<class SourceField, class TargetField>
void SharedDomainMap<CoordinateField>::apply( 
    const Teuchos::RCP< FieldEvaluator<MeshId,SourceField> >& source_evaluator,
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
    MeshId source_size = source_field_copy.size() / field_dim;
    Teuchos::RCP<Tpetra::MultiVector<typename SFT::value_type, int, MeshId> > 
	source_vector = Tpetra::createMultiVectorFromView( 
	    d_source_map, source_field_copy, source_size, field_dim );

    // Construct a view of the target space. Fill the target space with zeros
    // so that points we didn't map get some data.
    Teuchos::ArrayRCP<typename TFT::value_type> target_field_view(0,0);
    MeshId target_size = 0;
    if ( target_exists )
    {
	target_field_view = FieldTools<TargetField>::nonConstView( 
	    *target_space_manager->field() );

	target_size = target_field_view.size() / field_dim;

	FieldTools<TargetField>::putScalar( 
	    *target_space_manager->field(), 0.0 );

	DTK_CHECK( 
	    target_size == Teuchos::as<MeshId>(
		d_target_map->getNodeNumElements()) );
    }

    // Build a multivector for the target space.
    Teuchos::RCP<Tpetra::MultiVector<typename TFT::value_type, int, MeshId> > 
	target_vector =	Tpetra::createMultiVectorFromView( 
	    d_target_map, target_field_view, target_size, field_dim );

    // Move the data from the source decomposition to the target
    // decomposition.
    target_vector->doExport( *source_vector, *d_source_to_target_exporter, 
			     Tpetra::INSERT );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the target points that are in the rendezvous decomposition box.
 *
 * \param target_coord_manager The manager containing the target coordinates
 * to search the box with.
 *
 * \param target_ordinals The globally unique ordinals for the target
 * coordinates. 
 *
 * \param targets_in_box The global ordinals of the target coordinates in the
 * box. If a target point was not found in the box, return an invalid ordinal,
 * std::numeric_limits<MeshId>::max(), in its positition.
 */
template<class CoordinateField>
void SharedDomainMap<CoordinateField>::getTargetPointsInBox(
    const BoundingBox& box, const CoordinateField& target_coords,
    const Teuchos::ArrayView<const MeshId>& target_ordinals,
    Teuchos::Array<MeshId>& targets_in_box )
{
    Teuchos::ArrayRCP<const double> target_coords_view =
	FieldTools<CoordinateField>::view( target_coords );
    MeshId dim_size = 
	FieldTools<CoordinateField>::dimSize( target_coords );

    DTK_REQUIRE( dim_size == 
		 Teuchos::as<MeshId>(target_ordinals.size()) );

    targets_in_box.resize( dim_size );
    int field_dim = CFT::dim( target_coords );
    Teuchos::Array<double> target_point( field_dim );
    for ( MeshId n = 0; n < dim_size; ++n )
    {
	for ( int d = 0; d < field_dim; ++d )
	{
	    target_point[d] = target_coords_view[ dim_size*d + n ];
	}

	if ( box.pointInBox( target_point ) )
	{
	    targets_in_box[n] = target_ordinals[n];
	}
	else
	{
	    targets_in_box[n] = std::numeric_limits<MeshId>::max();
	}

	// If we're keeping track of the points not being mapped, add this
	// point's local index to the list if its not in the box.
	if ( d_store_missed_points && targets_in_box[n] == 
	     std::numeric_limits<MeshId>::max() )
	{
	    d_missed_points.push_back(n);
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_SHAREDDOMAINMAP_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_SharedDomainMap_def.hpp
//---------------------------------------------------------------------------//
