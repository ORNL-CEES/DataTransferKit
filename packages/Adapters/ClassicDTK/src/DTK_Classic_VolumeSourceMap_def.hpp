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
 * \file DTK_Classic_VolumeSourceMap_def.hpp
 * \author Stuart R. Slattery
 * \brief Volume source map definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_Classic_VOLUMESOURCEMAP_DEF_HPP
#define DTK_Classic_VOLUMESOURCEMAP_DEF_HPP

#include <algorithm>
#include <limits>
#include <set>

#include "DTK_DBC.hpp"
#include "DTK_FunctionSpace.hpp"
#include "DTK_FieldMultiVector.hpp"

#include "DTK_EntityCenteredShapeFunction.hpp"
#include "DTK_EntityCenteredField.hpp"
#include "DTK_Point.hpp"
#include "DTK_BasicGeometryLocalMap.hpp"

#include "DTK_ClassicGeometricEntity.hpp"
#include "DTK_ClassicGeometricEntityLocalMap.hpp"

#include "DTK_Classic_FieldTools.hpp"
#include "DTK_Classic_FieldTraits.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_Import.hpp>
#include <Tpetra_Distributor.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Export.hpp>

namespace DataTransferKit
{
namespace Classic
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param comm The communicator over which the map is generated.
 *
 * \param dimension The dimension of the map. This should be consistent with
 * all source and target objects (i.e. only 3 dimensional geometries will be
 * accepted with a 3 dimensional map). We need this here so we have a global
 * baseline for all objects that may or may not exist on all processes.
 * 
 * \param store_missed_points Set to true if it is desired to keep track of
 * the local target points missed during map generation. The default value is
 * false. 
 *
 * \param geometric_tolerance Tolerance used for point-in-geometry checks. The
 * default value is 1.0e-6.
 */
template<class Geometry, class GlobalOrdinal, class CoordinateField>
VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::VolumeSourceMap(
    const RCP_Comm& comm, const int dimension, 
    bool store_missed_points,
    const double geometric_tolerance )
    : d_comm( comm )
    , d_dimension( dimension )
    , d_store_missed_points( store_missed_points )
    , d_geometric_tolerance( geometric_tolerance )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Generate the volume source map.
 *
 * \param source_geometry_manager Source geometry in the volume source
 * problem. A null RCP is a valid argument. This will be the case when a
 * geometry manager is only constructed on a subset of the processes that the
 * shared domain map is constructed over. Note that the source geometry must
 * exist only on processes that reside within the VolumeSourceMap
 * communicator.
 *
 * \param target_coord_manager Target coordinates in the shared domain
 * problem. A null RCP is a valid argument. This will be the case when a field
 * manager is only constructed on a subset of the processes that the shared
 * domain map is constructed over. Note that the target coordinates must exist
 * only on processes that reside within the VolumeSourceMap communicator.
 */
template<class Geometry, class GlobalOrdinal, class CoordinateField>
void VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::setup( 
    const RCP_GeometryManager& source_geometry_manager, 
    const RCP_CoordFieldManager& target_coord_manager )
{
    // Create existence values for the managers.
    bool source_exists = true;
    if ( source_geometry_manager.is_null() ) source_exists = false;
    bool target_exists = true;
    if ( target_coord_manager.is_null() ) target_exists = false;

    // Create local to global process indexers for the managers.
    RCP_Comm source_comm;
    if ( source_exists )
    {
	source_comm = source_geometry_manager->comm();
    }
    RCP_Comm target_comm;
    if ( target_exists )
    {
	target_comm = target_coord_manager->comm();
    }
    d_source_indexer = DataTransferKit::CommIndexer( d_comm, source_comm );
    d_target_indexer = DataTransferKit::CommIndexer( d_comm, target_comm );

    // Check the source and target dimensions for consistency.
    if ( source_exists )
    {
	DTK_REQUIRE( source_geometry_manager->dim() == d_dimension );
    }

    if ( target_exists )
    {
	DTK_REQUIRE( CFT::dim( *target_coord_manager->field() ) 
			  == d_dimension );
    }
    
    // Build the domain space and map from the source information.
    // -----------------------------------------------------------

    // Create an entity set from the local source geometry. DTK entities are
    // stored as pointers to the classic geometry objects.
    d_source_entity_set =
	Teuchos::rcp( new DataTransferKit::BasicEntitySet(d_comm,d_dimension) );
    if ( source_exists )
    {
	// The classic API allowed for input source geometries to have
	// repeated global ids on different procs. The new API also allows for
	// this but we need to know which proc is the owning proc for parallel
	// uniqueness. Make globally unique ids even if sources are not
	// globally unique because we don't know with the version 1 API which
	// proc the user intended to be the owner. The map will handle the
	// averaging of multiple contributions.
	std::size_t local_num_source_geom =
	    source_geometry_manager->localNumGeometry();
	std::size_t global_num_source_geom =
	    source_geometry_manager->globalNumGeometry();
	Teuchos::RCP<const Tpetra::Map<int,DataTransferKit::EntityId> >
	    source_unique_id_map =
	    Tpetra::createContigMap<int,DataTransferKit::EntityId>(
		global_num_source_geom,
		local_num_source_geom,
		source_comm );
	Teuchos::ArrayView<const DataTransferKit::EntityId> unique_src_id_view =
	    source_unique_id_map->getNodeElementList();
	d_source_entity_ids.assign( unique_src_id_view.begin(),
				    unique_src_id_view.end() );
		
	// Create the entity set and extract the geometry centroids for future
	// evaluation.
	Teuchos::ArrayRCP<Geometry> source_geometry =
	    source_geometry_manager->geometry();
		Teuchos::ArrayRCP<GlobalOrdinal> source_gids =
		source_geometry_manager->gids();
	d_source_eval_ids.resize( local_num_source_geom );
	d_source_centroids.resize( local_num_source_geom*d_dimension );
	Teuchos::Array<double> source_centroid;
	for ( int i = 0; i < local_num_source_geom; ++i )
	{
	    d_source_entity_set->addEntity(
		DataTransferKit::ClassicGeometricEntity<Geometry>(
		    Teuchos::ptrFromRef(source_geometry[i]),
		    d_source_entity_ids[i],
		    d_comm->getRank())
		);
	    d_source_eval_ids[i] = source_gids[i];
	    source_centroid = GT::centroid( source_geometry[i] );
	    for ( int d = 0; d < d_dimension; ++d )
	    {
		d_source_centroids[d*local_num_source_geom + i] =
		    source_centroid[d];
	    }
	}
    }

    // Create a local map.
    Teuchos::RCP<DataTransferKit::ClassicGeometricEntityLocalMap<Geometry> >
	source_local_map = Teuchos::rcp(
	    new DataTransferKit::ClassicGeometricEntityLocalMap<Geometry>() );

    // Create a shape function.
    Teuchos::RCP<DataTransferKit::EntityCenteredShapeFunction>
	source_shape_function =	Teuchos::rcp(
	    new DataTransferKit::EntityCenteredShapeFunction() );

    // Create a function space for the source.
    DataTransferKit::FunctionSpace domain_function_space(
	d_source_entity_set, source_local_map, source_shape_function, Teuchos::null );

    // Make a map for the source vectors.
    Teuchos::RCP<const Tpetra::Map<int,DataTransferKit::SupportId> >
    	domain_vector_map = Tpetra::createNonContigMap<int,DataTransferKit::SupportId>(
	    d_source_entity_ids(), d_comm );
    
    // Build the target space and map from the target information.
    // -----------------------------------------------------------

    // Compute a unique global ordinal for each point in the coordinate field.
    Teuchos::Array<GlobalOrdinal> target_ordinals;
    computePointOrdinals( target_coord_manager, target_ordinals );

    // Create an entity set from the local target points.
    d_target_entity_set =
	Teuchos::rcp( new DataTransferKit::BasicEntitySet(d_comm,d_dimension) );
    if ( target_exists )
    {
	Teuchos::ArrayRCP<const typename CFT::value_type> coords_view =
	    FieldTools<CoordinateField>::view( *target_coord_manager->field() );
	Teuchos::Array<double> target_coords( d_dimension );
	int local_num_targets = target_ordinals.size();
	d_target_entity_ids.resize( local_num_targets );
	for ( int i = 0; i < local_num_targets; ++i )
	{
	    for ( int d = 0; d < d_dimension; ++d )
	    {
		target_coords[d] = coords_view[d*local_num_targets + i];
	    }
	    d_target_entity_set->addEntity(
		DataTransferKit::Point( target_ordinals[i],
					d_comm->getRank(),
					target_coords )
		);
	    d_target_entity_ids[i] = target_ordinals[i];
	}
    }

    // Create a local map.
    Teuchos::RCP<DataTransferKit::BasicGeometryLocalMap>
	target_local_map = Teuchos::rcp(
	    new DataTransferKit::BasicGeometryLocalMap() );

    // Create a shape function.
    Teuchos::RCP<DataTransferKit::EntityCenteredShapeFunction>
	target_shape_function =	Teuchos::rcp(
	    new DataTransferKit::EntityCenteredShapeFunction() );

    // Create a function space for the target.
    DataTransferKit::FunctionSpace range_function_space(
	d_target_entity_set, target_local_map, target_shape_function, Teuchos::null );

    // Make a map for the target vectors.
    Teuchos::RCP<const Tpetra::Map<int,DataTransferKit::SupportId> >
    	range_vector_map = Tpetra::createNonContigMap<int,DataTransferKit::SupportId>(
	    d_target_entity_ids(), d_comm );

    // Build the DTK operator.
    // -----------------------
    
    // Create parameters for the mapping.
    Teuchos::ParameterList parameters;
    Teuchos::ParameterList& search_list = parameters.sublist("Search");
    search_list.set<bool>("Track Missed Range Entities",d_store_missed_points);
    search_list.set<double>("Point Inclusion Tolerance", d_geometric_tolerance );
    
    // Create the interpolation operator.
    d_consistent_operator = Teuchos::rcp(
	new DataTransferKit::ConsistentInterpolationOperator<double>(
	    domain_vector_map, range_vector_map, parameters) );
    d_consistent_operator->setup( Teuchos::rcpFromRef(domain_function_space),
				  Teuchos::rcpFromRef(range_function_space) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply the volume source map for a valid source field evaluator and
 * target data space to the target points that were mapped.
 *
 * \param source_evaluator Function evaluator used to apply the mapping. This
 * FieldEvaluator must be valid for the source geometry used to generate the
 * map.
 *
 * \param target_space_manager Target space into which the function
 * evaluations will be written. Enough space must be allocated to hold
 * evaluations at all points in all dimensions of the field.
 */
template<class Geometry, class GlobalOrdinal, class CoordinateField>
template<class SourceField, class TargetField>
void VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::apply( 
    const Teuchos::RCP< FieldEvaluator<GlobalOrdinal,SourceField> >& source_evaluator,
    Teuchos::RCP< FieldManager<TargetField> >& target_space_manager )
{
    typedef FieldTraits<SourceField> SFT;
    typedef FieldTraits<TargetField> TFT;

    // Set existence values for the source and target.
    bool source_exists = true;
    if ( source_evaluator.is_null() ) source_exists = false;
    bool target_exists = true;
    if ( target_space_manager.is_null() ) target_exists = false;

    // Evaluate the source function at the centroids to extract the data.
    int source_dim;
    Teuchos::ArrayRCP<typename SFT::value_type> source_field_copy(0,0);
    if ( source_exists )
    {
	SourceField function_evaluations = 
	    source_evaluator->evaluate(
		Teuchos::arcpFromArray(d_source_eval_ids),
		Teuchos::arcpFromArray(d_source_centroids) );

	source_dim = SFT::dim( function_evaluations );

	source_field_copy =    
	    FieldTools<SourceField>::copy( function_evaluations );
    }
    Teuchos::broadcast<int,int>( *d_comm, d_source_indexer.l2g(0),
				 Teuchos::Ptr<int>(&source_dim) );

    // Build a vector for the source values.
    Teuchos::RCP<DataTransferKit::EntityCenteredField<typename SFT::value_type> >
	source_dtk_field = Teuchos::rcp(
	    new DataTransferKit::EntityCenteredField<typename SFT::value_type>(
		d_source_entity_ids(),
		source_dim,
		source_field_copy,
		DataTransferKit::EntityCenteredField<typename SFT::value_type>::BLOCKED)
	    );
    DataTransferKit::FieldMultiVector<typename SFT::value_type> source_vector(
	source_dtk_field, d_source_entity_set );

    // Construct a view of the target space.
    int target_dim;
    Teuchos::ArrayRCP<typename TFT::value_type> target_field_view(0,0);
    if ( target_exists )
    {
	target_field_view = FieldTools<TargetField>::nonConstView( 
	    *target_space_manager->field() );

	target_dim = TFT::dim( *target_space_manager->field() );
    }
    Teuchos::broadcast<int,int>( *d_comm, d_target_indexer.l2g(0),
				 Teuchos::Ptr<int>(&target_dim) );
    
    // Check that the source and target have the same field dimension.
    DTK_REQUIRE( source_dim == target_dim );

    // Verify that the target space has the proper amount of memory allocated.
    GlobalOrdinal target_size = target_field_view.size() / target_dim;
    if ( target_exists )
    {
	DTK_REQUIRE( target_size ==
		     Teuchos::as<GlobalOrdinal>(d_target_entity_ids.size()) );
    }

    // Build a vector for the target values.
    Teuchos::RCP<DataTransferKit::EntityCenteredField<typename TFT::value_type> >
	target_dtk_field = Teuchos::rcp(
	    new DataTransferKit::EntityCenteredField<typename TFT::value_type>(
		d_target_entity_ids(),
		target_dim,
		target_field_view,
		DataTransferKit::EntityCenteredField<typename TFT::value_type>::BLOCKED)
	    );
    DataTransferKit::FieldMultiVector<typename TFT::value_type> target_vector(
	target_dtk_field, d_target_entity_set );

    // Apply the DTK operator.
    d_consistent_operator->apply( source_vector, target_vector );
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
template<class Geometry, class GlobalOrdinal, class CoordinateField>
Teuchos::ArrayView<const GlobalOrdinal>
VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::getMissedTargetPoints() const
{
    DTK_REQUIRE( d_store_missed_points );
    Teuchos::ArrayView<const DataTransferKit::EntityId> missed_ids =
	d_consistent_operator->getMissedRangeEntityIds();
    d_missed_points.assign( missed_ids.begin(), missed_ids.end() );
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
template<class Geometry, class GlobalOrdinal, class CoordinateField>
Teuchos::ArrayView<GlobalOrdinal> 
VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::getMissedTargetPoints()
{
    DTK_REQUIRE( d_store_missed_points );
    Teuchos::ArrayView<const DataTransferKit::EntityId> missed_ids =
	d_consistent_operator->getMissedRangeEntityIds();
    d_missed_points.assign( missed_ids.begin(), missed_ids.end() );
    return d_missed_points();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute globally unique ordinals for the target points. Here an
 * invalid ordinal will be designated as the maximum value as specified by the
 * limits header for the ordinal type. We do this so that 0 may be a valid
 * ordinal.
 *
 * \param target_coords The coordinates to compute global ordinals for.
 *
 * \param target_ordinals The computed globally unique ordinals for the target
 * coordinates. 
 */
template<class Geometry, class GlobalOrdinal, class CoordinateField>
void 
VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::computePointOrdinals(
    const RCP_CoordFieldManager& target_coord_manager,
    Teuchos::Array<GlobalOrdinal>& target_ordinals )
{
    // Set an existence value for the target coords.
    bool target_exists = true;
    if ( target_coord_manager.is_null() ) target_exists = false;
    int comm_rank = d_comm->getRank();
    GlobalOrdinal local_size = 0;

    if ( target_exists )
    {
	int point_dim = CFT::dim( *target_coord_manager->field() );
	local_size = std::distance( 
	    CFT::begin( *target_coord_manager->field() ),
	    CFT::end( *target_coord_manager->field() ) ) / point_dim;
    }

    GlobalOrdinal global_max;
    Teuchos::reduceAll<int,GlobalOrdinal>( *d_comm,
					   Teuchos::REDUCE_MAX,
					   1,
					   &local_size,
					   &global_max );

    target_ordinals.resize( local_size );
    for ( GlobalOrdinal n = 0; n < local_size; ++n )
    {
	target_ordinals[n] = comm_rank*global_max + n;
    }
}

//---------------------------------------------------------------------------//

} // end namespace Classic
} // end namespace DataTransferKit

#endif // end DTK_Classic_VOLUMESOURCEMAP_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_Classic_VolumeSourceMap_def.hpp
//---------------------------------------------------------------------------//


