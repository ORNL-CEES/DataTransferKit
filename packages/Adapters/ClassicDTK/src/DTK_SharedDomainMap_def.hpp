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
#include <unordered_map>

#include "DTK_DBC.hpp"
#include "DTK_ParallelSearch.hpp"

#include "DTK_ClassicMesh.hpp"
#include "DTK_ClassicMeshEntitySet.hpp"
#include "DTK_ClassicMeshElementLocalMap.hpp"
#include "DTK_Point.hpp"
#include "DTK_BasicGeometryLocalMap.hpp"

#include "DTK_FieldTools.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Ptr.hpp>

#include <Tpetra_Import.hpp>
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
template<class Mesh, class CoordinateField>
SharedDomainMap<Mesh,CoordinateField>::SharedDomainMap(
    const RCP_Comm& comm, const int dimension, bool store_missed_points )
    : d_comm( comm )
    , d_dimension( dimension )
    , d_store_missed_points( store_missed_points )
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
 * \param tolerance Absolute tolerance for point searching. Will be used when
 * checking the reference cell ( and is therefore absolute ).
 */
template<class Mesh, class CoordinateField>
void SharedDomainMap<Mesh,CoordinateField>::setup(
    const RCP_MeshManager& source_mesh_manager,
    const RCP_CoordFieldManager& target_coord_manager,
    double tolerance )
{
    // Create existence values for the managers.
    bool source_exists = true;
    if ( source_mesh_manager.is_null() ) source_exists = false;
    bool target_exists = true;
    if ( target_coord_manager.is_null() ) target_exists = false;

    // Create local to global process indexers for the managers.
    RCP_Comm source_comm;
    if ( source_exists )
    {
        source_comm = source_mesh_manager->comm();
    }
    RCP_Comm target_comm;
    if ( target_exists )
    {
        target_comm = target_coord_manager->comm();
    }
    d_source_indexer = CommIndexer( d_comm, source_comm );
    d_target_indexer = CommIndexer( d_comm, target_comm );

    // Check the source and target dimensions for consistency.
    if ( source_exists )
    {
        DTK_REQUIRE( source_mesh_manager->dim() == d_dimension );
    }

    if ( target_exists )
    {
        DTK_REQUIRE( CFT::dim( *target_coord_manager->field() )
                          == d_dimension );
    }

    // Build the domain space and map from the source information.
    // -----------------------------------------------------------

    // Create an entity set from the local source mesh.
    Teuchos::RCP<DataTransferKit::ClassicMesh<Mesh> > classic_mesh =
        Teuchos::rcp( new DataTransferKit::ClassicMesh<Mesh>(source_mesh_manager) );
    ClassicMeshEntitySet<Mesh> source_entity_set( classic_mesh );

    // Create a local map.
    ClassicMeshElementLocalMap<Mesh> source_local_map(classic_mesh);

    // Build the target space and map from the target information.
    // -----------------------------------------------------------

    // Compute a unique global ordinal for each point in the coordinate field.
    Teuchos::Array<GlobalOrdinal> target_ordinals;
    computePointOrdinals( target_coord_manager, target_ordinals );

    // Create an entity set from the local target points.
    BasicEntitySet target_entity_set( d_comm, d_dimension );
    if ( target_exists )
    {
        Teuchos::ArrayRCP<const typename CFT::value_type> coords_view =
            FieldTools<CoordinateField>::view( *target_coord_manager->field() );
        Teuchos::Array<double> target_coords( d_dimension );
        int local_num_targets = target_ordinals.size();
        for ( int i = 0; i < local_num_targets; ++i )
        {
            for ( int d = 0; d < d_dimension; ++d )
            {
                target_coords[d] = coords_view[d*local_num_targets + i];
            }
            target_entity_set.addEntity(
                DataTransferKit::Point( target_ordinals[i],
                                        d_comm->getRank(),
                                        target_coords )
                );
        }
    }

    // Create a local map.
    DataTransferKit::BasicGeometryLocalMap target_local_map;

    // Find the location of the target points in the source mesh.
    // --------------------------------------------------------------

    // Create parameters for the mapping.
    Teuchos::ParameterList search_list;
    search_list.set<bool>("Track Missed Range Entities",d_store_missed_points);
    search_list.set<double>("Point Inclusion Tolerance", 1.0e-9 );

    // Do the parallel search.
    EntityIterator source_iterator = source_entity_set.entityIterator( d_dimension );
    EntityIterator target_iterator = target_entity_set.entityIterator( 0 );
    ParallelSearch parallel_search( d_comm,
                                    d_dimension,
                                    source_iterator,
                                    Teuchos::rcpFromRef(source_local_map),
                                    search_list );
    parallel_search.search( target_iterator,
                            Teuchos::rcpFromRef(target_local_map),
                            search_list );

    // Build the mapping.
    // -----------------------

    // Get the source-target parings.
    EntityIterator source_begin = source_iterator.begin();
    EntityIterator source_end = source_iterator.end();
    Teuchos::Array<EntityId> found_targets;
    Teuchos::Array<std::pair<EntityId,EntityId> > src_tgt_pairs;
    for ( auto src_geom = source_begin; src_geom != source_end; ++src_geom )
    {
        // Get the target points found in this source geometry.
        parallel_search.getRangeEntitiesFromDomain(
            src_geom->id(), found_targets );

        // If we found any points, add them to the mapping.
        for ( auto found_tgt : found_targets )
        {
            src_tgt_pairs.push_back(
                std::make_pair(src_geom->id(),found_tgt) );
        }
    }

    // Filter the source-target pairings so we only find a target point in one
    // geometry on this process. This handles the local uniqueness
    // problem. The tpetra import will handle the global uniqueness problem.
    auto sort_func = [] (std::pair<EntityId,EntityId> a,
                         std::pair<EntityId,EntityId> b )
                     { return a.second < b.second; };
    std::sort( src_tgt_pairs.begin(), src_tgt_pairs.end(), sort_func );
    auto unique_func = [] (std::pair<EntityId,EntityId> a,
                           std::pair<EntityId,EntityId> b )
                       { return a.second == b.second; };
    auto unique_it = std::unique( src_tgt_pairs.begin(),
                                  src_tgt_pairs.end(),
                                  unique_func );

    // Extract the mapping data.
    int num_tgt = std::distance( src_tgt_pairs.begin(), unique_it );
    Teuchos::Array<GlobalOrdinal> source_ordinals( num_tgt );
    d_source_geometry.resize( num_tgt );
    d_target_coords.resize( num_tgt * d_dimension );
    Teuchos::ArrayView<const double> tgt_coords;
    for ( int i = 0; i < num_tgt; ++i )
    {
        // Get the source geom id.
        d_source_geometry[i] = src_tgt_pairs[i].first;

        // Get the target point id.
        source_ordinals[i] = src_tgt_pairs[i].second;

        // Get the coordinates of the target point.
        parallel_search.rangeParametricCoordinatesInDomain(
            src_tgt_pairs[i].first,
            src_tgt_pairs[i].second,
            tgt_coords );

        for ( int d = 0; d < d_dimension; ++d )
        {
            d_target_coords[ d*num_tgt + i ] = tgt_coords[d];
        }
    }

    // Create the data map in the source decomposition.
    d_source_map = Tpetra::createNonContigMap<int,GlobalOrdinal>(
        source_ordinals(), d_comm );

    // Create the data map in the target decomposition.
    d_target_map = Tpetra::createNonContigMap<int,GlobalOrdinal>(
        target_ordinals(), d_comm );

    // Build the source-to-target importer.
    d_source_to_target_importer =
      Teuchos::rcp( new Tpetra::Import<int,GlobalOrdinal>(
          d_source_map, d_target_map ) );

    // Extract the missed points.
    if ( d_store_missed_points )
    {
        std::unordered_map<GlobalOrdinal,int> target_g2l;
        int local_num_targets = target_ordinals.size();
        for ( int t = 0; t < local_num_targets; ++t )
        {
            target_g2l.emplace( target_ordinals[t], t );
        }

        Teuchos::ArrayView<const EntityId> missed =
            parallel_search.getMissedRangeEntityIds();

        int num_missed = missed.size();
        d_missed_points.resize( num_missed );
        for ( int i = 0; i < num_missed; ++i )
        {
            DTK_CHECK( target_g2l.count(missed[i]) );
            d_missed_points[i] =
                target_g2l.find( missed[i] )->second;
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \Brief Apply the shared domain map for a valid source field evaluator and
 * target data space to the target points that were mapped.
 *
 * \param source_evaluator Function evaluator used to apply the mapping. This
 * FieldEvaluator must be valid for the source mesh used to generate the map.
 *
 * \param target_space_manager Target space into which the function
 * evaluations will be written. Enough space must be allocated to hold
 * evaluations at all points in all dimensions of the field.
 */
template<class Mesh, class CoordinateField>
template<class SourceField, class TargetField>
void SharedDomainMap<Mesh,CoordinateField>::apply(
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

    // Evaluate the source function at the target points and construct a view
    // of the function evaluations.
    int source_dim = 0;
    Teuchos::ArrayRCP<typename SFT::value_type> source_field_copy(0,0);
    if ( source_exists )
    {
        SourceField function_evaluations =
            source_evaluator->evaluate(
                Teuchos::arcpFromArray( d_source_geometry ),
                Teuchos::arcpFromArray( d_target_coords ) );

        source_dim = SFT::dim( function_evaluations );

        source_field_copy =
            FieldTools<SourceField>::copy( function_evaluations );
    }
    Teuchos::broadcast<int,int>( *d_comm, d_source_indexer.l2g(0),
                                 Teuchos::Ptr<int>(&source_dim) );

    // Build a multivector for the function evaluations.
    Tpetra::MultiVector<typename SFT::value_type, int, GlobalOrdinal>
        source_vector( d_source_map, source_dim );
    source_vector.get1dViewNonConst().deepCopy( source_field_copy() );

    // Construct a view of the target space.
    int target_dim = 0;
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
    if ( target_exists )
    {
        DTK_REQUIRE(
            target_field_view.size() == Teuchos::as<GlobalOrdinal>(
                d_target_map->getNodeNumElements()) * target_dim );
    }

    // Build a multivector for the target space.
    Tpetra::MultiVector<typename TFT::value_type, int, GlobalOrdinal>
        target_vector( d_target_map, target_dim );

    // Fill the target space with zeros so that points we didn't map get some
    // data.
    if ( target_exists )
    {
        FieldTools<TargetField>::putScalar(
            *target_space_manager->field(), 0.0 );
    }

    // Move the data from the source decomposition to the target
    // decomposition.
    target_vector.doImport( source_vector, *d_source_to_target_importer,
                            Tpetra::INSERT );
    target_field_view.deepCopy( target_vector.get1dView()() );
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
template<class Mesh, class CoordinateField>
Teuchos::ArrayView<const typename
                   SharedDomainMap<Mesh,CoordinateField>::GlobalOrdinal>
SharedDomainMap<Mesh,CoordinateField>::getMissedTargetPoints() const
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
template<class Mesh, class CoordinateField>
Teuchos::ArrayView<typename
                   SharedDomainMap<Mesh,CoordinateField>::GlobalOrdinal>
SharedDomainMap<Mesh,CoordinateField>::getMissedTargetPoints()
{
    DTK_REQUIRE( d_store_missed_points );
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
template<class Mesh, class CoordinateField>
void SharedDomainMap<Mesh,CoordinateField>::computePointOrdinals(
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

} // end namespace DataTransferKit

#endif // end DTK_SHAREDDOMAINMAP_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_SharedDomainMap_def.hpp
//---------------------------------------------------------------------------//
