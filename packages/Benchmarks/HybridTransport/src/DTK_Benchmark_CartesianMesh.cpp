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
/*!
 * \file DTK_Benchmark_CartesianMesh.hpp
 * \brief Local Cartesian mesh interface for the hybrid transport benchmark.
 */
//---------------------------------------------------------------------------//

#include "DTK_Benchmark_CartesianMesh.hpp"

namespace DataTransferKit
{
namespace Benchmark
{
//---------------------------------------------------------------------------//
// Constructor.
CartesianMesh::CartesianMesh(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm, const int set_id,
    const int block_id, const int num_i_blocks, const int num_j_blocks,
    const int num_k_blocks, const int x_global_num_node,
    const int y_global_num_node, const int x_edge_offset,
    const int y_edge_offset, const int z_edge_offset,
    const std::vector<double> &local_x_edges,
    const std::vector<double> &local_y_edges,
    const std::vector<double> &local_z_edges )
{
    // Set values.
    _comm = comm;
    _set_id = set_id;
    _block_id = block_id;
    _num_i_blocks = num_i_blocks;
    _num_j_blocks = num_j_blocks;
    _num_k_blocks = num_k_blocks;

    // Compute the local number of nodes.
    int x_local_num_node = local_x_edges.size();
    int y_local_num_node = local_y_edges.size();
    int z_local_num_node = local_z_edges.size();
    int local_num_node = x_local_num_node * y_local_num_node * z_local_num_node;

    // Local node indexer. Given a local ijk index compute a local cardinal
    // node index.
    auto local_node_index = [=]( const int i, const int j, const int k ) {
        return i + x_local_num_node * j +
               x_local_num_node * y_local_num_node * k;
    };

    // Global node indexer. Given a local ijk index compute a global cardinal
    // node index.
    auto global_node_index = [=]( const int i, const int j, const int k ) {
        int i_global = i + x_edge_offset;
        int j_global = j + y_edge_offset;
        int k_global = k + z_edge_offset;
        return i_global + x_global_num_node * j_global +
               x_global_num_node * y_global_num_node * k_global;
    };

    // Compute the local node global ids and coordinates.
    int space_dim = 3;
    _local_node_global_ids =
        Kokkos::View<GlobalOrdinal *>( "global_node_ids", local_num_node );
    auto local_node_global_ids_host =
        Kokkos::create_mirror_view( _local_node_global_ids );
    _local_node_coords =
        Kokkos::View<Coordinate **>( "node_coords", local_num_node, space_dim );
    auto local_node_coords_host =
        Kokkos::create_mirror_view( _local_node_coords );
    for ( int k = 0; k < z_local_num_node; ++k )
    {
        for ( int j = 0; j < y_local_num_node; ++j )
        {
            for ( int i = 0; i < x_local_num_node; ++i )
            {
                // Compute local and global indices.
                auto global_index = global_node_index( i, j, k );
                auto local_index = local_node_index( i, j, k );

                // Set the node global id.
                local_node_global_ids_host( local_index ) = global_index;

                // Get the node coordinates from the edge arrays.
                local_node_coords_host( local_index, 0 ) = local_x_edges[i];
                local_node_coords_host( local_index, 1 ) = local_y_edges[j];
                local_node_coords_host( local_index, 2 ) = local_z_edges[k];
            }
        }
    }
    Kokkos::deep_copy( _local_node_global_ids, local_node_global_ids_host );
    Kokkos::deep_copy( _local_node_coords, local_node_coords_host );

    // Compute the local number of cells.
    int x_local_num_cell = x_local_num_node - 1;
    int y_local_num_cell = y_local_num_node - 1;
    int z_local_num_cell = z_local_num_node - 1;
    int local_num_cell = x_local_num_cell * y_local_num_cell * z_local_num_cell;

    // Compute the global number of cells.
    int x_global_num_cell = x_global_num_node - 1;
    int y_global_num_cell = y_global_num_node - 1;

    // Local cell indexer. Given a local ijk index compute a local cardinal
    // cell index.
    auto local_cell_index = [=]( const int i, const int j, const int k ) {
        return i + x_local_num_cell * j +
               x_local_num_cell * y_local_num_cell * k;
    };

    // Global cell indexer. Given a local ijk index compute a global cardinal
    // cell index.
    auto global_cell_index = [=]( const int i, const int j, const int k ) {
        int i_global = i + x_edge_offset;
        int j_global = j + y_edge_offset;
        int k_global = k + z_edge_offset;
        return i_global + x_global_num_cell * j_global +
               x_global_num_cell * y_global_num_cell * k_global;
    };

    // Compute the local cell global ids, connectivities, and cell center
    // coordinates.
    int cell_num_node = 8;
    _local_cell_global_ids =
        Kokkos::View<GlobalOrdinal *>( "global_cell_ids", local_num_cell );
    auto local_cell_global_ids_host =
        Kokkos::create_mirror_view( _local_cell_global_ids );
    _local_cell_connectivity = Kokkos::View<LocalOrdinal **>(
        "cell_connectivity", local_num_cell, cell_num_node );
    auto local_cell_connectivity_host =
        Kokkos::create_mirror_view( _local_cell_connectivity );
    _local_cell_center_coords =
        Kokkos::View<Coordinate **>( "cell_coords", local_num_cell, space_dim );
    auto local_cell_center_coords_host =
        Kokkos::create_mirror_view( _local_cell_center_coords );
    for ( int k = 0; k < z_local_num_cell; ++k )
    {
        for ( int j = 0; j < y_local_num_cell; ++j )
        {
            for ( int i = 0; i < x_local_num_cell; ++i )
            {
                // Compute local and global indices.
                auto global_index = global_cell_index( i, j, k );
                auto local_index = local_cell_index( i, j, k );

                // Set the cell global id.
                local_cell_global_ids_host( local_index ) = global_index;

                // Set the cell connectivity. Connectivity is ordered for a
                // hex-8 cell in canonical ordering.
                local_cell_connectivity_host( local_index, 0 ) =
                    local_node_index( i, j, k );
                local_cell_connectivity_host( local_index, 1 ) =
                    local_node_index( i + 1, j, k );
                local_cell_connectivity_host( local_index, 2 ) =
                    local_node_index( i + 1, j + 1, k );
                local_cell_connectivity_host( local_index, 3 ) =
                    local_node_index( i, j + 1, k );
                local_cell_connectivity_host( local_index, 4 ) =
                    local_node_index( i, j, k + 1 );
                local_cell_connectivity_host( local_index, 5 ) =
                    local_node_index( i + 1, j, k + 1 );
                local_cell_connectivity_host( local_index, 6 ) =
                    local_node_index( i + 1, j + 1, k + 1 );
                local_cell_connectivity_host( local_index, 7 ) =
                    local_node_index( i, j + 1, k + 1 );

                // Set the cell center coordinates.
                local_cell_center_coords_host( local_index, 0 ) =
                    ( local_x_edges[i] + local_x_edges[i + 1] ) / 2.0;
                local_cell_center_coords_host( local_index, 1 ) =
                    ( local_y_edges[j] + local_y_edges[j + 1] ) / 2.0;
                local_cell_center_coords_host( local_index, 2 ) =
                    ( local_z_edges[k] + local_z_edges[k + 1] ) / 2.0;
            }
        }
    }
    Kokkos::deep_copy( _local_cell_global_ids, local_cell_global_ids_host );
    Kokkos::deep_copy( _local_cell_connectivity, local_cell_connectivity_host );
    Kokkos::deep_copy( _local_cell_center_coords,
                       local_cell_center_coords_host );
}

//---------------------------------------------------------------------------//

} // end namespace Benchmark
} // end namespace DataTransferKit
