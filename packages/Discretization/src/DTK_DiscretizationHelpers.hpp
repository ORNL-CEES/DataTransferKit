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

#ifndef DTK_DISCRETIZATION_HELPERS
#define DTK_DISCRETIZATION_HELPERS

#include <DTK_Topology.hpp>

#include <Kokkos_Macros.hpp>
#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
namespace Discretization
{
namespace Helpers
{
template <typename DeviceType>
std::array<unsigned int, DTK_N_TOPO> computeNCellsPerTopology(
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view )
{
    // This needs to be done on the host because we will use n_cells_per_topo
    // to allocate Kokkos::View and this needs to be done on the host
    std::array<unsigned int, DTK_N_TOPO> n_cells_per_topo;
    auto cell_topologies_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    Kokkos::deep_copy( cell_topologies_host, cell_topologies_view );
    unsigned int const n_cells = cell_topologies_view.extent( 0 );
    // Initialize dtk_cell_topo with an invalid topology
    DTK_CellTopology dtk_cell_topo = DTK_N_TOPO;
    n_cells_per_topo.fill( 0 );
    for ( unsigned int i = 0; i < n_cells; ++i )
    {
        dtk_cell_topo = cell_topologies_host( i );
        n_cells_per_topo[dtk_cell_topo] += 1;
    }

#if HAVE_DTK_DBC
    Topologies topologies;
    unsigned int dim = topologies[dtk_cell_topo].dim;
    // We do not support meshes that contain both 2D and 3D cells. All the
    // cells are either 2D or 3D
    for ( unsigned int i = 0; i < DTK_N_TOPO; ++i )
    {
        if ( n_cells_per_topo[i] != 0 )
            DTK_REQUIRE( topologies[dtk_cell_topo].dim == dim );
    }
#endif

    // PointInCell does not support all Intrepid2 topologies
    DTK_REQUIRE( n_cells_per_topo[DTK_TET_11] == 0 );
    DTK_REQUIRE( n_cells_per_topo[DTK_HEX_20] == 0 );
    DTK_REQUIRE( n_cells_per_topo[DTK_WEDGE_15] == 0 );

#if HAVE_DTK_DBC
    // The sum of the number of cells per topology should be equal to the size
    // of cell_topologies_view
    unsigned int sum = 0;
    for ( unsigned int i = 0; i < DTK_N_TOPO; ++i )
    {
        sum += n_cells_per_topo[i];
    }
    DTK_REQUIRE( sum == n_cells );
#endif

    return n_cells_per_topo;
}

template <typename DeviceType>
void checkOffsetOverflow( Kokkos::View<unsigned int *, DeviceType> offset )
{
#if HAVE_DTK_DBC
    // Check that we didn't overflow offset
    Kokkos::View<int[1], DeviceType> overflow( "overflow" );

    unsigned int const size = overflow.extent( 0 );
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::parallel_for(
        DTK_MARK_REGION( "check_overflow" ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, size - 1 ),
        KOKKOS_LAMBDA( int const i ) {
            if ( static_cast<int>( offset( i + 1 ) - offset( i ) ) < 0 )
                overflow[0] = 1;
        } );
    Kokkos::fence();
    auto overflow_host = Kokkos::create_mirror_view( overflow );
    Kokkos::deep_copy( overflow_host, overflow );
    DTK_REQUIRE( overflow_host( 0 ) == 0 );
#else
    (void)offset;
#endif
}

template <typename DeviceType, typename T1, typename T2>
void computeOffset( Kokkos::View<T1 *, DeviceType> predicate, T2 value,
                    Kokkos::View<unsigned int *, DeviceType> offset )
{
    DTK_REQUIRE( predicate.extent( 0 ) == offset.extent( 0 ) );

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

    ExecutionSpace space;
    ArborX::exclusivePrefixSum( space, mask, offset );

    // Check that we didn't overflow offset
    checkOffsetOverflow( offset );
}

template <typename DeviceType>
void computeNodeOffset(
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies,
    Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo,
    Kokkos::View<unsigned int *, DeviceType> node_offset )
{
    unsigned int const n_cells = cell_topologies.extent( 0 );
    Kokkos::View<unsigned int *, DeviceType> nodes_per_cell( "nodes_per_cell",
                                                             n_cells );

    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::parallel_for( DTK_MARK_REGION( "fill_nodes_per_cell" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_cells ),
                          KOKKOS_LAMBDA( int const i ) {
                              nodes_per_cell( i ) =
                                  n_nodes_per_topo( cell_topologies( i ) );
                          } );
    Kokkos::fence();

    ExecutionSpace space;
    ArborX::exclusivePrefixSum( space, nodes_per_cell, node_offset );

    // Check that we didn't overflow node_offset
    checkOffsetOverflow( node_offset );
}

template <typename DeviceType>
struct MeshOffsets
{
    MeshOffsets( Mesh<DeviceType> const &mesh )
        : n_nodes_per_topo( "n_nodes_per_topo" )
    {
        // First, we get the number of nodes for each topology to get
        // initialize block_cells at the right size.
        auto n_nodes_per_topo_host =
            Kokkos::create_mirror_view( n_nodes_per_topo );
        Topologies topologies;
        for ( int i = 0; i < DTK_N_TOPO; ++i )
            n_nodes_per_topo_host( i ) = topologies[i].n_nodes;
        Kokkos::deep_copy( n_nodes_per_topo, n_nodes_per_topo_host );

        unsigned int const n_cells = mesh.cell_topologies.extent( 0 );

        for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
        {
            offsets[topo_id] = Kokkos::View<unsigned int *, DeviceType>(
                "offset_" + std::to_string( topo_id ), n_cells );
            computeOffset( mesh.cell_topologies, topo_id, offsets[topo_id] );

            node_offsets[topo_id] = Kokkos::View<unsigned int *, DeviceType>(
                "node_offset", n_cells );
            computeNodeOffset( mesh.cell_topologies, n_nodes_per_topo,
                               node_offsets[topo_id] );
        }
    }

    std::array<Kokkos::View<unsigned int *, DeviceType>, DTK_N_TOPO> offsets;
    std::array<Kokkos::View<unsigned int *, DeviceType>, DTK_N_TOPO>
        node_offsets;
    Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo;
};

template <typename DeviceType>
KOKKOS_FUNCTION void
buildBlockCells( unsigned int const dim, int const i,
                 unsigned int const n_nodes, unsigned int const node_offset,
                 Kokkos::View<unsigned int *, DeviceType> cells,
                 Kokkos::View<unsigned int *, DeviceType> offset,
                 Kokkos::View<Coordinate **, DeviceType> coordinates,
                 Kokkos::View<Coordinate ***, DeviceType> block_cells )
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

/**
 * Convert the 1D Kokkos View cells and coordinates to arrays of 3D Kokkos
 * Views more suitable for Intrepid2.
 */
template <typename DeviceType>
void convertMesh( Mesh<DeviceType> const &mesh,
                  MeshOffsets<DeviceType> const &mesh_offsets,
                  std::array<Kokkos::View<Coordinate ***, DeviceType>,
                             DTK_N_TOPO> &block_cells )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        DTK_REQUIRE( mesh_offsets.offsets[topo_id].extent( 0 ) ==
                     mesh.cell_topologies.extent( 0 ) );
        DTK_REQUIRE( block_cells[topo_id].extent( 2 ) ==
                     mesh.nodes_coordinates.extent( 1 ) );

        unsigned int const dim = mesh.nodes_coordinates.extent( 1 );
        unsigned int const n_cells = mesh.cell_topologies.extent( 0 );
        auto block_cells_topo = block_cells[topo_id];
        auto node_offset = mesh_offsets.node_offsets[topo_id];
        auto offset = mesh_offsets.offsets[topo_id];

        Kokkos::parallel_for(
            DTK_MARK_REGION( "build_block_cells_" + std::to_string( topo_id ) ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_cells ),
            KOKKOS_LAMBDA( int const i ) {
                if ( mesh.cell_topologies( i ) == topo_id )
                {
                    buildBlockCells( dim, i,
                                     mesh_offsets.n_nodes_per_topo( topo_id ),
                                     node_offset( i ), mesh.cells, offset,
                                     mesh.nodes_coordinates, block_cells_topo );
                }
            } );
        Kokkos::fence();
    }
}

template <typename DeviceType>
KOKKOS_FUNCTION void
buildBoundingBoxes( unsigned int const dim, int const i,
                    unsigned int const n_nodes, unsigned int const node_offset,
                    Kokkos::View<unsigned int *, DeviceType> cells,
                    unsigned int const offset,
                    Kokkos::View<Coordinate **, DeviceType> coordinates,
                    Kokkos::View<Coordinate ***, DeviceType> block_cells,
                    Kokkos::View<ArborX::Box *, DeviceType> bounding_boxes )
{
    ArborX::Box bounding_box;
    // If dim == 2, we need to set bounding_box.minCorner()[2] and
    // bounding_box.maxCorner[2].
    if ( dim == 2 )
    {
        bounding_box.minCorner()[2] = 0;
        bounding_box.maxCorner()[2] = 1;
    }
    for ( unsigned int node = 0; node < n_nodes; ++node )
    {
        unsigned int const n = node_offset + node;
        for ( unsigned int d = 0; d < dim; ++d )
        {
            // Copy the coordinated in block_cells
            block_cells( offset, node, d ) = coordinates( cells( n ), d );
            // Build the bounding box.
            if ( block_cells( offset, node, d ) < bounding_box.minCorner()[d] )
                bounding_box.minCorner()[d] = block_cells( offset, node, d );
            if ( block_cells( offset, node, d ) > bounding_box.maxCorner()[d] )
                bounding_box.maxCorner()[d] = block_cells( offset, node, d );
        }
    }
    bounding_boxes( i ) = bounding_box;
}

/**
 * Build the bounding boxes associated to the cell and the map between the
 * bounding boxes and the flat array of cells.
 */
template <typename DeviceType>
void createBoundingBoxes(
    Mesh<DeviceType> const &mesh, MeshOffsets<DeviceType> const &mesh_offsets,
    std::array<Kokkos::View<Coordinate ***, DeviceType>, DTK_N_TOPO> const
        &block_cells,
    Kokkos::View<ArborX::Box *, DeviceType> bounding_boxes,
    Kokkos::View<unsigned int **, DeviceType> bounding_box_to_cell )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        DTK_REQUIRE( mesh_offsets.offsets[topo_id].extent( 0 ) ==
                     mesh.cell_topologies.extent( 0 ) );
        DTK_REQUIRE( mesh_offsets.node_offsets[topo_id].extent( 0 ) ==
                     mesh.cell_topologies.extent( 0 ) );
        DTK_REQUIRE( bounding_boxes.extent( 0 ) ==
                     mesh.cell_topologies.extent( 0 ) );

        // Build BoundingBoxes
        unsigned int const dim = mesh.nodes_coordinates.extent( 1 );
        unsigned int const n_cells = mesh.cell_topologies.extent( 0 );
        auto node_offset = mesh_offsets.node_offsets[topo_id];
        auto block_cells_topo = block_cells[topo_id];
        auto offset = mesh_offsets.offsets[topo_id];

        Kokkos::parallel_for(
            DTK_MARK_REGION( "build_bounding_boxes_" +
                             std::to_string( topo_id ) ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_cells ),
            KOKKOS_LAMBDA( int const i ) {
                if ( mesh.cell_topologies( i ) == topo_id )
                {
                    buildBoundingBoxes(
                        dim, i, mesh_offsets.n_nodes_per_topo( topo_id ),
                        node_offset( i ), mesh.cells, offset( i ),
                        mesh.nodes_coordinates, block_cells_topo,
                        bounding_boxes );
                }
            } );
        Kokkos::fence();

        // Build map between BoundingBoxes and BlockCells
        Kokkos::parallel_for(
            DTK_MARK_REGION( "build_bounding_boxes_to_block_cells_" +
                             std::to_string( topo_id ) ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_cells ),
            KOKKOS_LAMBDA( int const i ) {
                if ( mesh.cell_topologies( i ) == topo_id )
                {
                    bounding_box_to_cell( i, topo_id ) = offset( i );
                }
            } );
        Kokkos::fence();
    }
}
} // namespace Helpers
} // namespace Discretization
} // namespace DataTransferKit

#endif
