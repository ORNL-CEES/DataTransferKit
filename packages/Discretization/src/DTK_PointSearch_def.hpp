/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
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

#include <DTK_DBC.hpp>
#include <DTK_DetailsTeuchosSerializationTraits.hpp>
#include <DTK_DetailsUtils.hpp>
#include <DTK_DistributedSearchTree.hpp>
#include <DTK_PointInCell.hpp>
#include <DTK_Topology.hpp>

namespace DataTransferKit
{
namespace internal
{
template <typename DeviceType>
KOKKOS_FUNCTION void
buildBlockCells( unsigned int const dim, int const i,
                 unsigned int const n_nodes, unsigned int const node_offset,
                 Kokkos::View<unsigned int *, DeviceType> cells,
                 Kokkos::View<unsigned int *, DeviceType> offset,
                 Kokkos::View<double **, DeviceType> coordinates,
                 Kokkos::View<double ***, DeviceType> block_cells )
{
    unsigned int const k = offset( i );
    for ( unsigned int node = 0; node < n_nodes; ++node )
    {
        unsigned int const n = node_offset + node;
        for ( unsigned int d = 0; d < dim; ++d )
        {
            // Copy the coordinated in block_cells
            block_cells( k, node, d ) = coordinates( cells( n ), d );
        }
    }
}

template <typename DeviceType>
KOKKOS_FUNCTION void
buildBoundingBoxes( unsigned int const dim, int const i,
                    unsigned int const n_nodes, unsigned int const node_offset,
                    Kokkos::View<unsigned int *, DeviceType> cells,
                    Kokkos::View<unsigned int *, DeviceType> offset,
                    Kokkos::View<double **, DeviceType> coordinates,
                    Kokkos::View<double ***, DeviceType> block_cells,
                    Kokkos::View<Box *, DeviceType> bounding_boxes )
{
    Box bounding_box;
    // If dim == 2, we need to set bounding_box.minCorner()[2] and
    // bounding_box.maxCorner[2].
    if ( dim == 2 )
    {
        bounding_box.minCorner()[2] = 0;
        bounding_box.maxCorner()[2] = 1;
    }
    unsigned int const k = offset( i );
    for ( unsigned int node = 0; node < n_nodes; ++node )
    {
        unsigned int const n = node_offset + node;
        for ( unsigned int d = 0; d < dim; ++d )
        {
            // Copy the coordinated in block_cells
            block_cells( k, node, d ) = coordinates( cells( n ), d );
            // Build the bounding box.
            if ( block_cells( k, node, d ) < bounding_box.minCorner()[d] )
                bounding_box.minCorner()[d] = block_cells( k, node, d );
            if ( block_cells( k, node, d ) > bounding_box.maxCorner()[d] )
                bounding_box.maxCorner()[d] = block_cells( k, node, d );
        }
    }
    bounding_boxes( i ) = bounding_box;
}

template <typename DeviceType>
void convertPointDim( Kokkos::View<double **, DeviceType> points_coordinates,
                      unsigned int n_points,
                      Kokkos::View<double **, DeviceType> points_coord_3d )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::parallel_for( DTK_MARK_REGION( "convert_2D_pts_to_3D_pts" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_points ),
                          KOKKOS_LAMBDA( int const i ) {
                              points_coord_3d( i, 0 ) =
                                  points_coordinates( i, 0 );
                              points_coord_3d( i, 1 ) =
                                  points_coordinates( i, 1 );
                              points_coord_3d( i, 2 ) = 0.;
                          } );
    Kokkos::fence();
}

template <typename DeviceType>
void buildTopo( unsigned int const n_imports,
                Kokkos::View<int *, DeviceType> imported_cell_indices,
                Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell,
                Kokkos::View<unsigned int *, DeviceType> topo,
                Kokkos::View<unsigned int *, DeviceType> topo_size )
{
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
}

template <typename DeviceType, typename T1, typename T2>
void computeOffset( Kokkos::View<T1 *, DeviceType> predicate, T2 value,
                    Kokkos::View<unsigned int *, DeviceType> offset )
{
    // Create a Kokkos::View with ones where predicate matches value and
    // with zeros everywhere else.
    unsigned int const size = predicate.extent( 0 );
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::View<unsigned int *, DeviceType> mask( "mask", size );
    Kokkos::parallel_for( DTK_MARK_REGION( "compute_mask" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, size ),
                          KOKKOS_LAMBDA( int const i ) {
                              if ( predicate( i ) == value )
                                  mask( i ) = 1;
                              else
                                  mask( i ) = 0;
                          } );
    Kokkos::fence();

    exclusivePrefixSum( mask, offset );
}

template <typename DeviceType>
void computeNodeOffset(
    unsigned int const n_cells,
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
    Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo,
    Kokkos::View<unsigned int *, DeviceType> nodes_per_cell,
    Kokkos::View<unsigned int *, DeviceType> node_offset )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::parallel_for( DTK_MARK_REGION( "fill_nodes_per_cell" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_cells ),
                          KOKKOS_LAMBDA( int const i ) {
                              nodes_per_cell( i ) =
                                  n_nodes_per_topo( cell_topologies( i ) );
                          } );
    Kokkos::fence();

    exclusivePrefixSum( nodes_per_cell, node_offset );
}
} // namespace internal

template <typename DeviceType>
PointSearch<DeviceType>::PointSearch(
    Teuchos::RCP<const Teuchos::Comm<int>> comm,
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
    Kokkos::View<unsigned int *, DeviceType> cells,
    Kokkos::View<double **, DeviceType> cell_nodes_coordinates,
    Kokkos::View<double **, DeviceType> points_coordinates )
    : _comm( comm )
    , _target_to_source_distributor( _comm )
{
    // Initialize bounding_box_to_cell to an invalid state
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell(
        "bounding_box_to_cell", cell_topologies.extent( 0 ), DTK_N_TOPO );
    Kokkos::deep_copy( bounding_box_to_cell, static_cast<unsigned int>( -1 ) );

    // Compute the number of cells of each of the supported topologies.
    std::array<unsigned int, DTK_N_TOPO> n_cells_per_topo =
        computeNCellsPerTopology( cell_topologies );

    // Convert the cells and cell_nodes_coordinates View to block_cells
    std::array<Kokkos::View<double ***, DeviceType>, DTK_N_TOPO> block_cells;
    Kokkos::View<Box *, DeviceType> bounding_boxes(
        "bounding_boxes", cell_topologies.extent( 0 ) );
    convertMesh( n_cells_per_topo, cell_topologies, cells,
                 cell_nodes_coordinates, block_cells, bounding_boxes,
                 bounding_box_to_cell );

    // Perform the distributed search
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> per_topo_ranks;
    Kokkos::View<Point *, DeviceType> imported_points( "imported_points", 0 );
    Kokkos::View<int *, DeviceType> imported_query_ids( "imported_query_ids",
                                                        0 );
    Kokkos::View<int *, DeviceType> imported_cell_indices( "imported_indices",
                                                           0 );
    Kokkos::View<int *, DeviceType> ranks( "ranks", 0 );
    // The distributed search only supports dim == 3, so we need to convert the
    // 2D points to 3D points.
    if ( _dim == 2 )
    {
        unsigned int const n_points = points_coordinates.extent( 0 );
        Kokkos::View<double **, DeviceType> points_coord_3d( "points_coord_3d",
                                                             n_points, 3 );
        internal::convertPointDim<DeviceType>( points_coordinates, n_points,
                                               points_coord_3d );

        performDistributedSearch( points_coord_3d, bounding_boxes,
                                  imported_points, imported_query_ids,
                                  imported_cell_indices, ranks );
    }
    else
    {
        performDistributedSearch( points_coordinates, bounding_boxes,
                                  imported_points, imported_query_ids,
                                  imported_cell_indices, ranks );
    }

    // We need to separate the data for the different topologies because of
    // Intrepid2. Because a point can be found in multiple cells, we need to
    // compute the number of cells of each topology.
    unsigned int const n_imports = imported_points.extent( 0 );
    Kokkos::View<unsigned int *, DeviceType> topo( "topo", n_imports );
    Kokkos::View<unsigned int *, DeviceType> topo_size( "topo_size",
                                                        DTK_N_TOPO );
    internal::buildTopo( n_imports, imported_cell_indices, bounding_box_to_cell,
                         topo, topo_size );
    auto topo_size_host = Kokkos::create_mirror_view( topo_size );
    Kokkos::deep_copy( topo_size_host, topo_size );

    // Now that we know the size, allocate all the Views.
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO>
        filtered_per_topo_cell_indices;
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO>
        filtered_per_topo_query_ids;
    std::array<Kokkos::View<double **, DeviceType>, DTK_N_TOPO>
        filtered_per_topo_reference_points;
    std::array<Kokkos::View<bool *, DeviceType>, DTK_N_TOPO>
        filtered_per_topo_point_in_cell;
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO>
        filtered_per_topo_ranks;
    for ( unsigned int i = 0; i < DTK_N_TOPO; ++i )
    {
        filtered_per_topo_cell_indices[i] = Kokkos::View<int *, DeviceType>(
            "filtered_per_topo_cell_indices_" + std::to_string( i ),
            topo_size_host( i ) );
        filtered_per_topo_query_ids[i] = Kokkos::View<int *, DeviceType>(
            "filtered_per_topo_query_ids_" + std::to_string( i ),
            topo_size_host( i ) );
        filtered_per_topo_reference_points[i] =
            Kokkos::View<double **, DeviceType>(
                "filtered_per_topo_reference_points_" + std::to_string( i ),
                topo_size_host( i ), _dim );
        filtered_per_topo_point_in_cell[i] = Kokkos::View<bool *, DeviceType>(
            "filtered_per_topo_point_in_cell_" + std::to_string( i ),
            topo_size_host( i ) );
        filtered_per_topo_ranks[i] = Kokkos::View<int *, DeviceType>(
            "filtered_per_topo_ranks_" + std::to_string( i ),
            topo_size_host( i ) );
    }

    // Check if the points are in the cells
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
        if ( block_cells[topo_id].extent( 0 ) != 0 )
        {
            Kokkos::View<double **, DeviceType> filtered_points(
                "filtered_points", topo_size_host( topo_id ), _dim );
            performPointInCell( block_cells[topo_id], bounding_box_to_cell,
                                imported_cell_indices, imported_points,
                                imported_query_ids, ranks, topo, topo_id,
                                filtered_points,
                                filtered_per_topo_cell_indices[topo_id],
                                filtered_per_topo_query_ids[topo_id],
                                filtered_per_topo_reference_points[topo_id],
                                filtered_per_topo_point_in_cell[topo_id],
                                filtered_per_topo_ranks[topo_id] );
        }

    // Filter the points. Only keep the points that are in cell
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> filtered_ranks;
    filterInCell( filtered_per_topo_point_in_cell,
                  filtered_per_topo_reference_points,
                  filtered_per_topo_cell_indices, filtered_per_topo_query_ids,
                  filtered_per_topo_ranks, filtered_ranks );

    // Build the _source_to_target_distributor
    build_distributor( filtered_ranks );

    // Build a map between the cell_indices sorted by topology and the flat View
    // given to the constructor
    auto cell_topologies_host = Kokkos::create_mirror_view( cell_topologies );
    Kokkos::deep_copy( cell_topologies_host, cell_topologies );
    unsigned int const size = cell_topologies_host.extent( 0 );
    for ( unsigned int i = 0; i < size; ++i )
        _cell_indices_map[cell_topologies_host( i )].push_back( i );
}

template <typename DeviceType>
std::tuple<Kokkos::View<int *, DeviceType>, Kokkos::View<int *, DeviceType>,
           Kokkos::View<Point *, DeviceType>,
           Kokkos::View<unsigned int *, DeviceType>>
PointSearch<DeviceType>::getSearchResults()
{
    // Flatten the results
    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int n_ref_pts = 0;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
        n_ref_pts += _reference_points[topo_id].extent( 0 );

    Kokkos::View<int *, DeviceType> ranks( "ranks", n_ref_pts );
    Kokkos::deep_copy( ranks, _comm->getRank() );
    Kokkos::View<int *, DeviceType> cell_indices( "cell_indices", n_ref_pts );
    auto cell_indices_host = Kokkos::create_mirror_view( cell_indices );
    Kokkos::View<unsigned int *, DeviceType> query_ids( "query_ids",
                                                        n_ref_pts );
    Kokkos::View<Point *, DeviceType> ref_pts( "ref_pts", n_ref_pts );
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
                                      ref_pts( i + n_copied_pts )[d] =
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
    Kokkos::View<Point *, DeviceType> imported_ref_pts( "imported_ref_pts",
                                                        n_imports );
    Kokkos::View<unsigned int *, DeviceType> imported_query_ids(
        "imported_query_ids", n_imports );

    Details::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        _target_to_source_distributor, ranks, imported_ranks );
    Details::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        _target_to_source_distributor, cell_indices, imported_cell_indices );
    Details::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        _target_to_source_distributor, ref_pts, imported_ref_pts );
    Details::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        _target_to_source_distributor, query_ids, imported_query_ids );

    Details::DistributedSearchTreeImpl<DeviceType>::sortResults(
        imported_query_ids, imported_query_ids, imported_cell_indices,
        imported_ranks, imported_ref_pts );

    return std::make_tuple( imported_ranks, imported_cell_indices,
                            imported_ref_pts, imported_query_ids );
}

template <typename DeviceType>
void PointSearch<DeviceType>::buildBlockCells(
    unsigned int n_cells, unsigned int topo_id,
    std::array<Kokkos::View<double ***, DeviceType>, DTK_N_TOPO> const
        &block_cells,
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
    Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo,
    Kokkos::View<unsigned int *, DeviceType> node_offset,
    Kokkos::View<unsigned int *, DeviceType> cells,
    Kokkos::View<unsigned int *, DeviceType> offset,
    Kokkos::View<double **, DeviceType> coordinates )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int const dim = _dim;
    Kokkos::View<double ***, DeviceType> block_cells_topo =
        block_cells[topo_id];
    Kokkos::parallel_for(
        DTK_MARK_REGION( "build_block_cells_" + std::to_string( topo_id ) ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, n_cells ),
        KOKKOS_LAMBDA( int const i ) {
            if ( cell_topologies( i ) == topo_id )
            {
                internal::buildBlockCells( dim, i, n_nodes_per_topo( topo_id ),
                                           node_offset( i ), cells, offset,
                                           coordinates, block_cells_topo );
            }
        } );
    Kokkos::fence();
}

template <typename DeviceType>
void PointSearch<DeviceType>::buildBoundingBoxes(
    unsigned int n_cells, unsigned int topo_id,
    std::array<Kokkos::View<double ***, DeviceType>, DTK_N_TOPO> const
        &block_cells,
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
    Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo,
    Kokkos::View<unsigned int *, DeviceType> node_offset,
    Kokkos::View<unsigned int *, DeviceType> cells,
    Kokkos::View<unsigned int *, DeviceType> offset,
    Kokkos::View<double **, DeviceType> coordinates,
    Kokkos::View<Box *, DeviceType> bounding_boxes )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int const dim = _dim;
    Kokkos::View<double ***, DeviceType> block_cells_topo =
        block_cells[topo_id];
    Kokkos::parallel_for(
        DTK_MARK_REGION( "build_bounding_boxes_" + std::to_string( topo_id ) ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, n_cells ),
        KOKKOS_LAMBDA( int const i ) {
            if ( cell_topologies( i ) == topo_id )
            {
                internal::buildBoundingBoxes(
                    dim, i, n_nodes_per_topo( topo_id ), node_offset( i ),
                    cells, offset, coordinates, block_cells_topo,
                    bounding_boxes );
            }
        } );
    Kokkos::fence();
}

template <typename DeviceType>
void PointSearch<DeviceType>::buildBoundingBoxesToBlockCells(
    unsigned int n_cells, unsigned int topo_id,
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
    Kokkos::View<unsigned int *, DeviceType> offset,
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::parallel_for(
        DTK_MARK_REGION( "build_bounding_boxes_to_block_cells_" +
                         std::to_string( topo_id ) ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, n_cells ),
        KOKKOS_LAMBDA( int const i ) {
            if ( cell_topologies( i ) == topo_id )
            {
                bounding_box_to_cell( i, topo_id ) = offset( i );
            }
        } );
    Kokkos::fence();
}

template <typename DeviceType>
void PointSearch<DeviceType>::convertMesh(
    std::array<unsigned int, DTK_N_TOPO> const &n_cells_per_topo,
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
    Kokkos::View<unsigned int *, DeviceType> cells,
    Kokkos::View<double **, DeviceType> coordinates,
    std::array<Kokkos::View<double ***, DeviceType>, DTK_N_TOPO> &block_cells,
    Kokkos::View<Box *, DeviceType> bounding_boxes,
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell )
{
    // First, we get the number of nodes for each topology to get
    // initialize block_cells at the right size.
    Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo(
        "n_nodes_per_topo" );
    auto n_nodes_per_topo_host = Kokkos::create_mirror_view( n_nodes_per_topo );
    Topologies topologies;
    for ( int i = 0; i < DTK_N_TOPO; ++i )
        n_nodes_per_topo_host( i ) = topologies[i].n_nodes;
    Kokkos::deep_copy( n_nodes_per_topo, n_nodes_per_topo_host );

    // Create an array of Kokkos::View where each View is a block of
    // cells with the same topology
    for ( int i = 0; i < DTK_N_TOPO; ++i )
    {
        block_cells[i] = Kokkos::View<double ***, DeviceType>(
            "block_cells_" + std::to_string( i ), n_cells_per_topo[i],
            n_nodes_per_topo_host( i ), _dim );
    }

    // Compute the offset associated to each cell in the coordinates
    // View.
    unsigned int const n_cells = cell_topologies.extent( 0 );
    Kokkos::View<unsigned int *, DeviceType> nodes_per_cell( "nodes_per_cell",
                                                             n_cells );
    Kokkos::View<unsigned int *, DeviceType> node_offset( "node_offset",
                                                          n_cells );
    internal::computeNodeOffset( n_cells, cell_topologies, n_nodes_per_topo,
                                 nodes_per_cell, node_offset );

    Kokkos::View<unsigned int *, DeviceType> offset( "offset", n_cells );
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        // Compute the relative position of each cell of a given
        // topology
        internal::computeOffset( cell_topologies, topo_id, offset );

        // Build BlockCells
        buildBlockCells( n_cells, topo_id, block_cells, cell_topologies,
                         n_nodes_per_topo, node_offset, cells, offset,
                         coordinates );

        // Build BoundingBoxes
        buildBoundingBoxes( n_cells, topo_id, block_cells, cell_topologies,
                            n_nodes_per_topo, node_offset, cells, offset,
                            coordinates, bounding_boxes );

        // Build map between BoundingBoxes and BlockCells
        buildBoundingBoxesToBlockCells( n_cells, topo_id, cell_topologies,
                                        offset, bounding_box_to_cell );
    }
}

template <typename DeviceType>
void PointSearch<DeviceType>::performDistributedSearch(
    Kokkos::View<double **, DeviceType> points_coord,
    Kokkos::View<Box *, DeviceType> bounding_boxes,
    Kokkos::View<Point *, DeviceType> &imported_points,
    Kokkos::View<int *, DeviceType> &imported_query_ids,
    Kokkos::View<int *, DeviceType> &imported_cell_indices,
    Kokkos::View<int *, DeviceType> &ranks )
{
    DistributedSearchTree<DeviceType> distributed_tree( _comm, bounding_boxes );

    unsigned int const n_points = points_coord.extent( 0 );

    // Build the queries
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::View<Within *, DeviceType> queries( "queries", n_points );
    Kokkos::parallel_for( DTK_MARK_REGION( "register_queries" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_points ),
                          KOKKOS_LAMBDA( int i ) {
                              queries( i ) = within(
                                  {{points_coord( i, 0 ), points_coord( i, 1 ),
                                    points_coord( i, 2 )}},
                                  0. );
                          } );
    Kokkos::fence();

    // Perform the distributed search
    Kokkos::View<int *, DeviceType> indices( "indices" );
    Kokkos::View<int *, DeviceType> offset( "offset" );
    distributed_tree.query( queries, indices, offset, ranks );

    // Create the source to target distributor
    auto ranks_host = Kokkos::create_mirror_view( ranks );
    Kokkos::deep_copy( ranks_host, ranks );
    Tpetra::Distributor source_to_target_distributor( _comm );
    unsigned int const n_imports = source_to_target_distributor.createFromSends(
        Teuchos::ArrayView<int const>( ranks_host.data(),
                                       ranks_host.extent( 0 ) ) );

    // Communicate cell indices
    Kokkos::realloc( imported_cell_indices, n_imports );
    Details::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        source_to_target_distributor, indices, imported_cell_indices );
    // Duplicate the points_coord for the communication. Duplicating the points
    // allows us to use the same distributor.
    unsigned int const indices_size = indices.extent( 0 );
    Kokkos::View<Point *, DeviceType> exported_points( "exported_points",
                                                       indices_size );
    Kokkos::View<int *, DeviceType> exported_query_ids( "exported_query_ids",
                                                        indices_size );
    // This line should not be necessary but there is problem with the
    // lambda capture on CUDA otherwise.
    unsigned int dim = _dim;
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

    // Communicate the points
    Kokkos::realloc( imported_points, n_imports );
    Details::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        source_to_target_distributor, exported_points, imported_points );

    // Communicate the query_ids. We communicate the query_ids to keep track of
    // which points is associated to which query.
    Kokkos::realloc( imported_query_ids, n_imports );
    Details::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        source_to_target_distributor, exported_query_ids, imported_query_ids );

    // Communicate the ranks of the sending processors. This will be used to
    // build the _target_to_source_distributor.
    Kokkos::realloc( ranks, n_imports );
    Kokkos::View<int *, DeviceType> exported_ranks( "exported_ranks",
                                                    indices_size );
    Kokkos::deep_copy( exported_ranks, _comm->getRank() );
    Details::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        source_to_target_distributor, exported_ranks, ranks );
}

template <typename DeviceType>
void PointSearch<DeviceType>::filterTopology(
    Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell,
    Kokkos::View<int *, DeviceType> cell_indices,
    Kokkos::View<Point *, DeviceType> points,
    Kokkos::View<int *, DeviceType> query_ids,
    Kokkos::View<int *, DeviceType> ranks,
    Kokkos::View<int *, DeviceType> filtered_cell_indices,
    Kokkos::View<double **, DeviceType> filtered_points,
    Kokkos::View<int *, DeviceType> filtered_query_ids,
    Kokkos::View<int *, DeviceType> filtered_ranks )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int const n_imports = topo.extent( 0 );
    Kokkos::View<unsigned int *, DeviceType> offset( "offset", n_imports );
    internal::computeOffset( topo, topo_id, offset );

    // Create Kokkos::View with the points and the cell indices associated
    // with cells of topo_id topology. Also transform 3D points back to 2D
    // points.
    unsigned int dim = _dim;
    Kokkos::parallel_for(
        "filter_data", Kokkos::RangePolicy<ExecutionSpace>( 0, n_imports ),
        KOKKOS_LAMBDA( int const i ) {
            if ( topo( i ) == topo_id )
            {
                unsigned int const k = offset( i );
                filtered_cell_indices( k ) =
                    bounding_box_to_cell( cell_indices( i ), topo_id );
                for ( unsigned int j = 0; j < dim; ++j )
                    filtered_points( k, j ) = points( i )[j];
                filtered_query_ids( k ) = query_ids( i );
                filtered_ranks( k ) = ranks( i );
            }
        } );
    Kokkos::fence();
}

template <typename DeviceType>
void PointSearch<DeviceType>::filterInCell(
    std::array<Kokkos::View<bool *, DeviceType>, DTK_N_TOPO> const
        &filtered_per_topo_point_in_cell,
    std::array<Kokkos::View<double **, DeviceType>, DTK_N_TOPO> const
        &filtered_per_topo_reference_points,
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const
        &filtered_per_topo_cell_indices,
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const
        &filtered_per_topo_query_ids,
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> const &ranks,
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> &filtered_ranks )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int dim = _dim;

    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        unsigned int n_ref_points =
            filtered_per_topo_point_in_cell[topo_id].extent( 0 );
        if ( n_ref_points != 0 )
        {
            int n_filtered_ref_points = 0;
            Kokkos::View<bool *, DeviceType> pt_in_cell =
                filtered_per_topo_point_in_cell[topo_id];
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
            Kokkos::realloc( filtered_ranks[topo_id], n_filtered_ref_points );

            // We cannot use private member in a lambda function with CUDA
            Kokkos::View<Coordinate **, DeviceType> ref_points =
                _reference_points[topo_id];
            Kokkos::View<int *, DeviceType> query_ids = _query_ids[topo_id];
            Kokkos::View<int *, DeviceType> cell_indices =
                _cell_indices[topo_id];

            Kokkos::View<unsigned int *, DeviceType> offset( "offset",
                                                             n_ref_points );
            internal::computeOffset( pt_in_cell, true, offset );
            Kokkos::parallel_for(
                DTK_MARK_REGION( "filter" ),
                Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_points ),
                KOKKOS_LAMBDA( int const i ) {
                    if ( pt_in_cell[i] )
                    {
                        unsigned int k = offset( i );
                        for ( unsigned int d = 0; d < dim; ++d )
                            ref_points( k, d ) =
                                filtered_per_topo_reference_points[topo_id](
                                    i, d );
                        query_ids( k ) =
                            filtered_per_topo_query_ids[topo_id]( i );
                        cell_indices( k ) =
                            filtered_per_topo_cell_indices[topo_id]( i );
                        filtered_ranks[topo_id]( k ) = ranks[topo_id]( i );
                    }
                } );
            Kokkos::fence();
        }
    }
}

template <typename DeviceType>
std::array<unsigned int, DTK_N_TOPO>
PointSearch<DeviceType>::computeNCellsPerTopology(
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view )
{
    // This needs to be done on the host because we will use n_cells_per_topo
    // to allocate Kokkos::View and this needs to be done on the host
    std::array<unsigned int, DTK_N_TOPO> n_cells_per_topo;
    auto cell_topologies_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    Kokkos::deep_copy( cell_topologies_host, cell_topologies_view );
    unsigned int const n_cells = cell_topologies_view.extent( 0 );
    DTK_CellTopology dtk_cell_topo;
    n_cells_per_topo.fill( 0 );
    for ( unsigned int i = 0; i < n_cells; ++i )
    {
        dtk_cell_topo = cell_topologies_host( i );
        n_cells_per_topo[dtk_cell_topo] += 1;
    }

    // We do not support meshes that contains both 2D and 3D cells. All the
    // cells are either 2D or 3D
    Topologies topologies;
    _dim = topologies[dtk_cell_topo].dim;
#if HAVE_DTK_DBC
    for ( unsigned int i = 0; i < DTK_N_TOPO; ++i )
    {
        if ( n_cells_per_topo[i] != 0 )
            DTK_REQUIRE( topologies[dtk_cell_topo].dim == _dim );
    }
#endif

    // PointInCell does not support all Intrepid2 topologies
    DTK_REQUIRE( n_cells_per_topo[DTK_TET_11] == 0 );
    DTK_REQUIRE( n_cells_per_topo[DTK_HEX_20] == 0 );
    DTK_REQUIRE( n_cells_per_topo[DTK_WEDGE_15] == 0 );

    return n_cells_per_topo;
}

template <typename DeviceType>
void PointSearch<DeviceType>::performPointInCell(
    Kokkos::View<double ***, DeviceType> cells,
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell,
    Kokkos::View<int *, DeviceType> imported_cell_indices,
    Kokkos::View<Point *, DeviceType> imported_points,
    Kokkos::View<int *, DeviceType> imported_query_ids,
    Kokkos::View<int *, DeviceType> imported_ranks,
    Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
    Kokkos::View<double **, DeviceType> filtered_points,
    Kokkos::View<int *, DeviceType> filtered_cell_indices,
    Kokkos::View<int *, DeviceType> filtered_query_ids,
    Kokkos::View<double **, DeviceType> reference_points,
    Kokkos::View<bool *, DeviceType> point_in_cell,
    Kokkos::View<int *, DeviceType> filtered_ranks )
{
    filterTopology( topo, topo_id, bounding_box_to_cell, imported_cell_indices,
                    imported_points, imported_query_ids, imported_ranks,
                    filtered_cell_indices, filtered_points, filtered_query_ids,
                    filtered_ranks );

    // Perform the PointInCell search
    Topologies topologies;
    PointInCell<DeviceType>::search(
        filtered_points, cells, filtered_cell_indices, topologies[topo_id].topo,
        reference_points, point_in_cell );
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

    _target_to_source_distributor.createFromSends(
        Teuchos::ArrayView<int const>( flatten_ranks ) );
}
} // namespace DataTransferKit

// Explicit instantiation macro
#define DTK_POINTSEARCH_INSTANT( NODE )                                        \
    template class PointSearch<typename NODE::device_type>;

#endif
