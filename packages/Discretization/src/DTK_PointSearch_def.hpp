/****************************************************************************
 * Copyright (c) 2012-2020 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef DTK_POINT_SEARCH_DEF_HPP
#define DTK_POINT_SEARCH_DEF_HPP

#include <ArborX.hpp>
#include <DTK_DBC.hpp>
#include <DTK_DetailsUtils.hpp>
#include <DTK_DiscretizationHelpers.hpp>
#include <DTK_PointInCell.hpp>
#include <DTK_Topology.hpp>

#include <mpi.h>

namespace DataTransferKit
{
namespace internal
{
template <typename DeviceType>
Kokkos::View<Coordinate **, DeviceType>
convertPointDim( Kokkos::View<Coordinate **, DeviceType> points_coord_2d )
{
    DTK_REQUIRE( points_coord_2d.extent( 1 ) == 2 );

    unsigned int const n_points = points_coord_2d.extent( 0 );
    Kokkos::View<Coordinate **, DeviceType> points_coord_3d( "points_coord_3d",
                                                             n_points, 3 );

    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::parallel_for( DTK_MARK_REGION( "convert_2D_pts_to_3D_pts" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_points ),
                          KOKKOS_LAMBDA( int const i ) {
                              points_coord_3d( i, 0 ) = points_coord_2d( i, 0 );
                              points_coord_3d( i, 1 ) = points_coord_2d( i, 1 );
                              points_coord_3d( i, 2 ) = 0.;
                          } );
    Kokkos::fence();

    return points_coord_3d;
}

template <typename DeviceType>
void buildTopo( Kokkos::View<int *, DeviceType> imported_cell_indices,
                Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell,
                Kokkos::View<unsigned int *, DeviceType> topo,
                Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> topo_size )
{
    DTK_REQUIRE( bounding_box_to_cell.extent( 1 ) == DTK_N_TOPO );
    DTK_REQUIRE( topo.extent( 0 ) == imported_cell_indices.extent( 0 ) );

    unsigned int const n_imports = imported_cell_indices.extent( 0 );
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::parallel_for(
        DTK_MARK_REGION( "build_topo" ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, n_imports ),
        KOKKOS_LAMBDA( int const i ) {
            for ( unsigned int j = 0; j < DTK_N_TOPO; ++j )
            {
                if ( bounding_box_to_cell( imported_cell_indices( i ), j ) !=
                     static_cast<unsigned int>( -1 ) )
                {
                    topo( i ) = j;
                    Kokkos::atomic_increment( &topo_size( j ) );
                    break;
                }
            }
        } );
    Kokkos::fence();

#if HAVE_DTK_DBC
    // The sum of the values in topo_size should be equal to the size of
    // imported_cell_indices
    Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> topo_size_sum(
        "topo_size_sum" );
    ArborX::exclusivePrefixSum( topo_size, topo_size_sum );
    DTK_REQUIRE( ArborX::lastElement( topo_size_sum ) == n_imports );
#endif
}

template <typename ViewType>
void sendDataAcrossNetwork(
    ArborX::Details::Distributor<typename ViewType::device_type> const
        &distributor,
    std::pair<ViewType, ViewType> data )
{
    ArborX::Details::DistributedTreeImpl<typename ViewType::device_type>::
        sendAcrossNetwork( typename ViewType::execution_space{}, distributor,
                           data.first, data.second );
}

template <typename T, typename... Targs>
void sendDataAcrossNetwork(
    ArborX::Details::Distributor<typename T::device_type> const &distributor,
    std::pair<T, T> d, Targs... data )
{
    sendDataAcrossNetwork( distributor, d );
    sendDataAcrossNetwork( distributor, data... );
}

//  Return parameters points, cell_indices, query_ids,
template <typename DeviceType>
std::tuple<Kokkos::View<ArborX::Point *, DeviceType>,
           Kokkos::View<int *, DeviceType>, Kokkos::View<int *, DeviceType>,
           Kokkos::View<int *, DeviceType>>
moveDataFromSourceToTarget(
    MPI_Comm comm, Kokkos::View<int *, DeviceType> indices,
    Kokkos::View<int *, DeviceType> offset,
    Kokkos::View<int *, DeviceType> ranks,
    Kokkos::View<Coordinate **, DeviceType> points_coord, unsigned int dim )
{
    using ExecutionSpace = typename DeviceType::execution_space;

    // Create the source to target distributor
    ArborX::Details::Distributor<DeviceType> source_to_target_distributor(
        comm );
    unsigned int const n_imports =
        source_to_target_distributor.createFromSends( ExecutionSpace{}, ranks );

    // Duplicate the points_coord for the communication. Duplicating the points
    // allows us to use the same distributor.
    unsigned int const indices_size = indices.extent( 0 );
    Kokkos::View<ArborX::Point *, DeviceType> exported_points(
        "exported_points", indices_size );
    Kokkos::View<int *, DeviceType> exported_query_ids( "exported_query_ids",
                                                        indices_size );
    Kokkos::parallel_for(
        "duplicate_points",
        Kokkos::RangePolicy<ExecutionSpace>( 0, offset.extent( 0 ) - 1 ),
        KOKKOS_LAMBDA( int const i ) {
            for ( int j = offset( i ); j < offset( i + 1 ); ++j )
            {
                exported_query_ids( j ) = i;
                for ( unsigned int k = 0; k < dim; ++k )
                    exported_points( j )[k] = points_coord( i, k );
            }
        } );
    Kokkos::fence();

    Kokkos::View<int *, DeviceType> exported_ranks( "exported_ranks",
                                                    indices_size );
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    Kokkos::deep_copy( exported_ranks, comm_rank );

    Kokkos::View<ArborX::Point *, DeviceType> imported_points(
        "imported_points", n_imports );
    Kokkos::View<int *, DeviceType> imported_cell_indices( "imported_indices",
                                                           n_imports );
    Kokkos::View<int *, DeviceType> imported_query_ids( "imported_query_ids",
                                                        n_imports );
    Kokkos::View<int *, DeviceType> imported_ranks( "ranks", n_imports );

    sendDataAcrossNetwork(
        source_to_target_distributor,
        std::make_pair( exported_points, imported_points ),
        std::make_pair( indices, imported_cell_indices ),
        std::make_pair( exported_query_ids, imported_query_ids ),
        std::make_pair( exported_ranks, imported_ranks ) );

    return std::make_tuple( imported_points, imported_cell_indices,
                            imported_query_ids, imported_ranks );
}
} // namespace internal

template <typename DeviceType>
PointSearch<DeviceType>::PointSearch(
    MPI_Comm comm, Mesh<DeviceType> const &mesh,
    Kokkos::View<Coordinate **, DeviceType> points_coordinates )
    : _comm( comm )
    , _target_to_source_distributor( _comm )
{
    DTK_REQUIRE( points_coordinates.extent( 1 ) ==
                 mesh.nodes_coordinates.extent( 1 ) );
    _dim = points_coordinates.extent( 1 );

    // Compute the number of cells of each of the supported topologies.
    std::array<unsigned int, DTK_N_TOPO> n_cells_per_topo =
        Discretization::Helpers::computeNCellsPerTopology(
            mesh.cell_topologies );

    // Compute the topology and node offset
    Discretization::Helpers::MeshOffsets<DeviceType> mesh_offsets( mesh );

    // Convert the cells and cell_nodes_coordinates View to block_cells
    auto n_nodes_per_topo_host =
        Kokkos::create_mirror_view( mesh_offsets.n_nodes_per_topo );
    Kokkos::deep_copy( n_nodes_per_topo_host, mesh_offsets.n_nodes_per_topo );
    std::array<Kokkos::View<Coordinate ***, DeviceType>, DTK_N_TOPO>
        block_cells;
    for ( int i = 0; i < DTK_N_TOPO; ++i )
    {
        block_cells[i] = Kokkos::View<Coordinate ***, DeviceType>(
            "block_cells_" + std::to_string( i ), n_cells_per_topo[i],
            n_nodes_per_topo_host( i ), _dim );
    }
    Discretization::Helpers::convertMesh( mesh, mesh_offsets, block_cells );

    // Initialize bounding_box_to_cell to an invalid state
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell(
        "bounding_box_to_cell", mesh.cell_topologies.extent( 0 ), DTK_N_TOPO );
    Kokkos::deep_copy( bounding_box_to_cell, static_cast<unsigned int>( -1 ) );

    Kokkos::View<ArborX::Box *, DeviceType> bounding_boxes(
        "bounding_boxes", mesh.cell_topologies.extent( 0 ) );
    Discretization::Helpers::createBoundingBoxes(
        mesh, mesh_offsets, block_cells, bounding_boxes, bounding_box_to_cell );

    // Perform the distributed search. At the end of the distributed search the
    // points are moved from the "source processors" to the "target processors".
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> per_topo_ranks;
    Kokkos::View<ArborX::Point *, DeviceType> imported_points;
    Kokkos::View<int *, DeviceType> imported_query_ids;
    Kokkos::View<int *, DeviceType> imported_cell_indices;
    Kokkos::View<int *, DeviceType> imported_ranks;
    std::tie( imported_points, imported_cell_indices, imported_query_ids,
              imported_ranks ) =
        performDistributedSearch(
            ( _dim == 3 ) ? points_coordinates
                          : internal::convertPointDim( points_coordinates ),
            bounding_boxes );

    // We need to separate the data for the different topologies because of
    // Intrepid2. Because a point can be found in multiple cells, we need to
    // compute the number of cells of each topology.
    unsigned int const n_imports = imported_points.extent( 0 );
    Kokkos::View<unsigned int *, DeviceType> topo( "topo", n_imports );
    Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> topo_size( "topo_size" );
    internal::buildTopo( imported_cell_indices, bounding_box_to_cell, topo,
                         topo_size );
    auto topo_size_host = Kokkos::create_mirror_view( topo_size );
    Kokkos::deep_copy( topo_size_host, topo_size );

    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> filtered_ranks;
    // Check if the points are in the cells
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
        if ( block_cells[topo_id].extent( 0 ) != 0 )
        {
            filtered_ranks[topo_id] = performPointInCell(
                block_cells[topo_id], bounding_box_to_cell,
                imported_cell_indices, imported_points, imported_query_ids,
                imported_ranks, topo, topo_id, topo_size_host( topo_id ) );
        }

    // Build the _source_to_target_distributor
    build_distributor( filtered_ranks );

    // Build a map between the cell_indices sorted by topology and the flat View
    // given to the constructor
    auto cell_topologies_host =
        Kokkos::create_mirror_view( mesh.cell_topologies );
    Kokkos::deep_copy( cell_topologies_host, mesh.cell_topologies );
    unsigned int const size = cell_topologies_host.extent( 0 );
    for ( unsigned int i = 0; i < size; ++i )
        _cell_indices_map[cell_topologies_host( i )].push_back( i );
}

template <typename DeviceType>
std::tuple<Kokkos::View<int *, DeviceType>, Kokkos::View<int *, DeviceType>,
           Kokkos::View<Coordinate *[3], DeviceType>,
           Kokkos::View<unsigned int *, DeviceType>>
PointSearch<DeviceType>::getSearchResults() const {
    // Flatten the results
    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int n_ref_pts = 0;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
        n_ref_pts += _reference_points[topo_id].extent( 0 );

    Kokkos::View<int *, DeviceType> ranks( "ranks", n_ref_pts );
    int comm_rank;
    MPI_Comm_rank( _comm, &comm_rank );
    Kokkos::deep_copy( ranks, comm_rank );
    Kokkos::View<int *, DeviceType> cell_indices( "cell_indices", n_ref_pts );
    auto cell_indices_host = Kokkos::create_mirror_view( cell_indices );
    Kokkos::View<unsigned int *, DeviceType> query_ids( "query_ids",
                                                        n_ref_pts );
    Kokkos::View<Coordinate * [3], DeviceType> ref_pts( "ref_pts", n_ref_pts );
    unsigned int n_copied_pts = 0;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        unsigned int const size = _query_ids[topo_id].extent( 0 );

        // First fill cell_indices. This has to be done on the host because
        // _cell_indices_map only exists on the host.
        auto topo_cell_indices_host =
            Kokkos::create_mirror_view( _cell_indices[topo_id] );
        Kokkos::deep_copy( topo_cell_indices_host, _cell_indices[topo_id] );
        for ( unsigned int i = 0; i < size; ++i )
        {
            cell_indices_host( i + n_copied_pts ) =
                _cell_indices_map[topo_id][topo_cell_indices_host( i )];
        }

        // Fill query_ids
        auto topo_query_ids = _query_ids[topo_id];
        Kokkos::parallel_for( DTK_MARK_REGION( "query_ids" ),
                              Kokkos::RangePolicy<ExecutionSpace>( 0, size ),
                              KOKKOS_LAMBDA( int const i ) {
                                  query_ids( i + n_copied_pts ) =
                                      topo_query_ids( i );
                              } );
        Kokkos::fence();

        // Fill ref_pts
        unsigned int dim = _dim;
        auto topo_ref_pts = _reference_points[topo_id];
        Kokkos::parallel_for( DTK_MARK_REGION( "fill_ref_pts" ),
                              Kokkos::RangePolicy<ExecutionSpace>( 0, size ),
                              KOKKOS_LAMBDA( int const i ) {
                                  for ( unsigned int d = 0; d < dim; ++d )
                                      ref_pts( i + n_copied_pts, d ) =
                                          topo_ref_pts( i, d );
                              } );
        Kokkos::fence();

        n_copied_pts += size;
    }
    Kokkos::deep_copy( cell_indices, cell_indices_host );

    // Communicate the results
    unsigned int n_imports =
        _target_to_source_distributor.getTotalReceiveLength();
    Kokkos::View<int *, DeviceType> imported_ranks( "imported_ranks",
                                                    n_imports );
    Kokkos::View<int *, DeviceType> imported_cell_indices(
        "imported_cell_indices", n_imports );
    Kokkos::View<Coordinate * [3], DeviceType> imported_ref_pts(
        "imported_ref_pts", n_imports );
    Kokkos::View<unsigned int *, DeviceType> imported_query_ids(
        "imported_query_ids", n_imports );

    internal::sendDataAcrossNetwork(
        _target_to_source_distributor, std::make_pair( ranks, imported_ranks ),
        std::make_pair( cell_indices, imported_cell_indices ),
        std::make_pair( ref_pts, imported_ref_pts ),
        std::make_pair( query_ids, imported_query_ids ) );

    ArborX::Details::DistributedTreeImpl<DeviceType>::sortResults(
        ExecutionSpace{}, imported_query_ids, imported_query_ids,
        imported_cell_indices, imported_ranks, imported_ref_pts );

#if HAVE_DTK_DBC
    // Check that ranks and cell indices are positive
    Kokkos::View<int[2], DeviceType> negative_values( "negative_values" );

    Kokkos::parallel_for(
        DTK_MARK_REGION( "check_positivity" ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, n_imports ),
        KOKKOS_LAMBDA( int const i ) {
            if ( imported_ranks( i ) < 0 )
                Kokkos::atomic_increment( &negative_values( 0 ) );
            if ( imported_cell_indices( i ) < 0 )
                Kokkos::atomic_increment( &negative_values( 1 ) );
        } );
    Kokkos::fence();
    auto negative_values_host = Kokkos::create_mirror_view( negative_values );
    Kokkos::deep_copy( negative_values_host, negative_values );
    DTK_REQUIRE( negative_values_host( 0 ) == 0 );
    DTK_REQUIRE( negative_values_host( 1 ) == 0 );
#endif

    return std::make_tuple( imported_ranks, imported_cell_indices,
                            imported_ref_pts, imported_query_ids );
}

template <typename DeviceType>
std::tuple<Kokkos::View<ArborX::Point *, DeviceType>,
           Kokkos::View<int *, DeviceType>, Kokkos::View<int *, DeviceType>,
           Kokkos::View<int *, DeviceType>> PointSearch<DeviceType>::
    performDistributedSearch(
        Kokkos::View<Coordinate **, DeviceType> points_coord,
        Kokkos::View<ArborX::Box *, DeviceType> bounding_boxes )
{
    DTK_REQUIRE( points_coord.extent( 1 ) == 3 );

    using ExecutionSpace = typename DeviceType::execution_space;
    using MemorySpace = typename DeviceType::memory_space;

    ArborX::DistributedTree<MemorySpace> distributed_tree(
        _comm, ExecutionSpace{}, bounding_boxes );

    unsigned int const n_points = points_coord.extent( 0 );

    // Build the queries
    Kokkos::View<decltype( ArborX::intersects( ArborX::Sphere{} ) ) *,
                 DeviceType>
        queries( "queries", n_points );
    Kokkos::parallel_for( DTK_MARK_REGION( "register_queries" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_points ),
                          KOKKOS_LAMBDA( int i ) {
                              queries( i ) = ArborX::intersects( ArborX::Sphere{
                                  {static_cast<float>( points_coord( i, 0 ) ),
                                   static_cast<float>( points_coord( i, 1 ) ),
                                   static_cast<float>( points_coord( i, 2 ) )},
                                  0.} );
                          } );
    Kokkos::fence();

    // Perform the distributed search
    using PairIndexRank = Kokkos::pair<int, int>;
    Kokkos::View<int *, DeviceType> offset( "offset", 0 );
    Kokkos::View<PairIndexRank *, DeviceType> index_rank( "index_rank", 0 );
    distributed_tree.query( ExecutionSpace{}, queries, index_rank, offset );

    // Split the pair
    Kokkos::View<int *, DeviceType> indices( "indices", 0 );
    Kokkos::View<int *, DeviceType> ranks( "ranks", 0 );
    Details::splitIndexRank( index_rank, indices, ranks );

    // Move the points from the source processors to the target processors
    return internal::moveDataFromSourceToTarget( _comm, indices, offset, ranks,
                                                 points_coord, _dim );
}

template <typename DeviceType>
std::tuple<Kokkos::View<int *, DeviceType>,
           Kokkos::View<Coordinate **, DeviceType>,
           Kokkos::View<int *, DeviceType>, Kokkos::View<int *, DeviceType>>
PointSearch<DeviceType>::filterTopology(
    Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
    unsigned int size,
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell,
    Kokkos::View<int *, DeviceType> cell_indices,
    Kokkos::View<ArborX::Point *, DeviceType> points,
    Kokkos::View<int *, DeviceType> query_ids,
    Kokkos::View<int *, DeviceType> ranks )
{
    DTK_REQUIRE( topo.extent( 0 ) == ranks.extent( 0 ) );
    DTK_REQUIRE( query_ids.extent( 0 ) == ranks.extent( 0 ) );
    DTK_REQUIRE( bounding_box_to_cell.extent( 1 ) > topo_id );

    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int const n_imports = topo.extent( 0 );
    Kokkos::View<unsigned int *, DeviceType> offset( "offset", n_imports );
    Discretization::Helpers::computeOffset( topo, topo_id, offset );

    // Create Kokkos::View with the points and the cell indices associated
    // with cells of topo_id topology. Also transform 3D points back to 2D
    // points.
    Kokkos::View<Coordinate **, DeviceType> filtered_per_topo_points(
        "filtered_per_topo_points", size, _dim );
    Kokkos::View<int *, DeviceType> filtered_per_topo_cell_indices(
        "filtered_per_topo_cell_indices_" + std::to_string( topo_id ), size );
    Kokkos::View<int *, DeviceType> filtered_per_topo_query_ids(
        "filtered_per_topo_query_ids_" + std::to_string( topo_id ), size );
    Kokkos::View<int *, DeviceType> filtered_per_topo_ranks(
        "filtered_per_topo_ranks_" + std::to_string( topo_id ), size );
    unsigned int dim = _dim;
    Kokkos::parallel_for(
        "filter_data", Kokkos::RangePolicy<ExecutionSpace>( 0, n_imports ),
        KOKKOS_LAMBDA( int const i ) {
            if ( topo( i ) == topo_id )
            {
                unsigned int const k = offset( i );
                filtered_per_topo_cell_indices( k ) =
                    bounding_box_to_cell( cell_indices( i ), topo_id );
                for ( unsigned int j = 0; j < dim; ++j )
                    filtered_per_topo_points( k, j ) = points( i )[j];
                filtered_per_topo_query_ids( k ) = query_ids( i );
                filtered_per_topo_ranks( k ) = ranks( i );
            }
        } );
    Kokkos::fence();

    return std::make_tuple(
        filtered_per_topo_cell_indices, filtered_per_topo_points,
        filtered_per_topo_query_ids, filtered_per_topo_ranks );
}

template <typename DeviceType>
Kokkos::View<int *, DeviceType> PointSearch<DeviceType>::filterInCell(
    Kokkos::View<bool *, DeviceType> filtered_per_topo_point_in_cell,
    Kokkos::View<Coordinate **, DeviceType> filtered_per_topo_reference_points,
    Kokkos::View<int *, DeviceType> filtered_per_topo_cell_indices,
    Kokkos::View<int *, DeviceType> filtered_per_topo_query_ids,
    Kokkos::View<int *, DeviceType> filtered_per_topo_ranks,
    unsigned int topo_id )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int dim = _dim;

    Kokkos::View<int *, DeviceType> filtered_ranks;
    unsigned int n_ref_points = filtered_per_topo_point_in_cell.extent( 0 );
    if ( n_ref_points != 0 )
    {
        int n_filtered_ref_points = 0;
        Kokkos::View<bool *, DeviceType> pt_in_cell =
            filtered_per_topo_point_in_cell;
        Kokkos::parallel_reduce(
            DTK_MARK_REGION( "compute_n_ref_pts" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_points ),
            KOKKOS_LAMBDA( int i, int &partial_sum ) {
                if ( pt_in_cell[i] == true )
                    partial_sum += 1;
            },
            n_filtered_ref_points );

        // We are only interested in points that belong to the cells. So we
        // need to filter out all the points that were false positive of
        // the distributed search.
        Kokkos::realloc( _reference_points[topo_id], n_filtered_ref_points,
                         _dim );
        Kokkos::realloc( _query_ids[topo_id], n_filtered_ref_points );
        Kokkos::realloc( _cell_indices[topo_id], n_filtered_ref_points );
        Kokkos::realloc( filtered_ranks, n_filtered_ref_points );

        // We cannot use private member in a lambda function with CUDA
        Kokkos::View<Coordinate **, DeviceType> ref_points =
            _reference_points[topo_id];
        Kokkos::View<int *, DeviceType> query_ids = _query_ids[topo_id];
        Kokkos::View<int *, DeviceType> cell_indices = _cell_indices[topo_id];

        Kokkos::View<unsigned int *, DeviceType> offset( "offset",
                                                         n_ref_points );
        Discretization::Helpers::computeOffset( pt_in_cell, true, offset );
        Kokkos::parallel_for(
            DTK_MARK_REGION( "filter" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_points ),
            KOKKOS_LAMBDA( int const i ) {
                if ( pt_in_cell[i] )
                {
                    unsigned int k = offset( i );
                    for ( unsigned int d = 0; d < dim; ++d )
                        ref_points( k, d ) =
                            filtered_per_topo_reference_points( i, d );
                    query_ids( k ) = filtered_per_topo_query_ids( i );
                    cell_indices( k ) = filtered_per_topo_cell_indices( i );
                    filtered_ranks( k ) = filtered_per_topo_ranks( i );
                }
            } );
        Kokkos::fence();
    }

    return filtered_ranks;
}

template <typename DeviceType>
Kokkos::View<int *, DeviceType> PointSearch<DeviceType>::performPointInCell(
    Kokkos::View<Coordinate ***, DeviceType> cells,
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell,
    Kokkos::View<int *, DeviceType> imported_cell_indices,
    Kokkos::View<ArborX::Point *, DeviceType> imported_points,
    Kokkos::View<int *, DeviceType> imported_query_ids,
    Kokkos::View<int *, DeviceType> imported_ranks,
    Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
    unsigned int size )
{
    // Filter the data for a given topology
    Kokkos::View<Coordinate **, DeviceType> filtered_per_topo_points;
    Kokkos::View<int *, DeviceType> filtered_per_topo_cell_indices;
    Kokkos::View<int *, DeviceType> filtered_per_topo_query_ids;
    Kokkos::View<int *, DeviceType> filtered_per_topo_ranks;
    std::tie( filtered_per_topo_cell_indices, filtered_per_topo_points,
              filtered_per_topo_query_ids, filtered_per_topo_ranks ) =
        filterTopology( topo, topo_id, size, bounding_box_to_cell,
                        imported_cell_indices, imported_points,
                        imported_query_ids, imported_ranks );

    // Perform the PointInCell search
    Topologies topologies;
    Kokkos::View<Coordinate **, DeviceType> filtered_per_topo_reference_points(
        "filtered_per_topo_reference_points_" + std::to_string( topo_id ), size,
        _dim );
    Kokkos::View<bool *, DeviceType> filtered_per_topo_point_in_cell(
        "filtered_per_topo_point_in_cell_" + std::to_string( topo_id ), size );
    PointInCell<DeviceType>::search(
        filtered_per_topo_points, cells, filtered_per_topo_cell_indices,
        topologies[topo_id].topo, filtered_per_topo_reference_points,
        filtered_per_topo_point_in_cell );

    // Filter the points. Only keep the points that are in cell
    Kokkos::View<int *, DeviceType> filtered_ranks;
    filtered_ranks = filterInCell(
        filtered_per_topo_point_in_cell, filtered_per_topo_reference_points,
        filtered_per_topo_cell_indices, filtered_per_topo_query_ids,
        filtered_per_topo_ranks, topo_id );

    return filtered_ranks;
}

template <typename DeviceType>
void PointSearch<DeviceType>::build_distributor(
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const
        &filtered_ranks )
{
    // Flatten the filtered ranks to be used by the distributor
    std::vector<int> flatten_ranks;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        auto rank_host = Kokkos::create_mirror_view( filtered_ranks[topo_id] );
        Kokkos::deep_copy( rank_host, filtered_ranks[topo_id] );
        unsigned int const rank_host_size = rank_host.size();
        for ( unsigned int i = 0; i < rank_host_size; ++i )
            flatten_ranks.push_back( rank_host( i ) );
    }

    Kokkos::DefaultHostExecutionSpace space;
    _target_to_source_distributor.createFromSends(
        space, Kokkos::View<int const *, Kokkos::HostSpace,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                   flatten_ranks.data(), flatten_ranks.size() ) );
}
} // namespace DataTransferKit

// Explicit instantiation macro
#define DTK_POINTSEARCH_INSTANT( NODE )                                        \
    template class PointSearch<typename NODE::device_type>;

#endif
