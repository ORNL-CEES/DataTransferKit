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
 * \file DTK_Benchmark_MonteCarloMesh.cpp
 * \brief MonteCarlo mesh interface for the hybrid transport benchmark.
 */
//---------------------------------------------------------------------------//

#include "DTK_Benchmark_MonteCarloMesh.hpp"
#include "DTK_DBC.hpp"

#include <algorithm>

namespace DataTransferKit
{
namespace Benchmark
{
//---------------------------------------------------------------------------//
// Uniform cell size constructor.
MonteCarloMesh::MonteCarloMesh(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm, const int num_sets,
    const int num_cells_i, const int num_cells_j, const int num_cells_k,
    const double delta_x, const double delta_y, const double delta_z,
    const std::vector<double> &x_bnd_mesh,
    const std::vector<double> &y_bnd_mesh,
    const std::vector<double> &z_bnd_mesh )
{
    // Create uniform global edge arrays.
    std::vector<double> global_x_edges( num_cells_i + 1 );
    for ( int n = 0; n < num_cells_i + 1; ++n )
        global_x_edges[n] = n * delta_x;
    std::vector<double> global_y_edges( num_cells_j + 1 );
    for ( int n = 0; n < num_cells_j + 1; ++n )
        global_y_edges[n] = n * delta_y;
    std::vector<double> global_z_edges( num_cells_k + 1 );
    for ( int n = 0; n < num_cells_k + 1; ++n )
        global_z_edges[n] = n * delta_z;

    // Partition the mesh.
    partition( comm, num_sets, global_x_edges, global_y_edges, global_z_edges,
               x_bnd_mesh, y_bnd_mesh, z_bnd_mesh );
}

//---------------------------------------------------------------------------//
// Global edge constructor.
MonteCarloMesh::MonteCarloMesh(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm, const int num_sets,
    const std::vector<double> &global_x_edges,
    const std::vector<double> &global_y_edges,
    const std::vector<double> &global_z_edges,
    const std::vector<double> &x_bnd_mesh,
    const std::vector<double> &y_bnd_mesh,
    const std::vector<double> &z_bnd_mesh )
{
    // Partition the mesh.
    partition( comm, num_sets, global_x_edges, global_y_edges, global_z_edges,
               x_bnd_mesh, y_bnd_mesh, z_bnd_mesh );
}

//---------------------------------------------------------------------------//
// Partition the mesh.
void MonteCarloMesh::partition(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm, const int num_sets,
    const std::vector<double> &global_x_edges,
    const std::vector<double> &global_y_edges,
    const std::vector<double> &global_z_edges,
    const std::vector<double> &x_bnd_mesh,
    const std::vector<double> &y_bnd_mesh,
    const std::vector<double> &z_bnd_mesh )
{
    // Determine how many I, J, and K blocks there will be.
    int num_i_blocks = x_bnd_mesh.size() - 1;
    int num_j_blocks = y_bnd_mesh.size() - 1;
    int num_k_blocks = z_bnd_mesh.size() - 1;
    int num_blocks = num_i_blocks * num_j_blocks * num_k_blocks;
    DTK_REQUIRE( num_sets * num_blocks == comm->getSize() );

    // Calculate my set id and my block ids. All blocks in a single set are
    // given a contiguous group of ranks.
    int comm_rank = comm->getRank();

    int set_id = std::floor( comm_rank / num_blocks );
    DTK_CHECK( 0 <= set_id && set_id <= num_sets );

    int block_id = comm_rank % num_blocks;
    DTK_CHECK( 0 <= block_id && block_id <= num_blocks );

    // Calculate cardinal block ids.
    int k_block = std::floor( block_id / ( num_i_blocks * num_j_blocks ) );
    DTK_CHECK( 0 <= k_block && k_block < num_k_blocks );

    int j_block = std::floor(
        ( block_id - k_block * num_i_blocks * num_j_blocks ) / num_i_blocks );
    DTK_CHECK( 0 <= j_block && j_block < num_j_blocks );

    int i_block = block_id - j_block * num_i_blocks -
                  k_block * num_i_blocks * num_j_blocks;
    DTK_CHECK( 0 <= i_block && i_block < num_i_blocks );

    DTK_CHECK( block_id == i_block + j_block * num_i_blocks +
                               k_block * num_i_blocks * num_j_blocks );

    // Calculate the local edges based on the boundary mesh locations.
    std::vector<double> local_x_edges;
    int x_offset;
    computeLocalEdges( global_x_edges, x_bnd_mesh, i_block, local_x_edges,
                       x_offset );

    std::vector<double> local_y_edges;
    int y_offset;
    computeLocalEdges( global_y_edges, y_bnd_mesh, j_block, local_y_edges,
                       y_offset );

    std::vector<double> local_z_edges;
    int z_offset;
    computeLocalEdges( global_z_edges, z_bnd_mesh, k_block, local_z_edges,
                       z_offset );

    // Build the mesh data.
    _cartesian_mesh = std::make_shared<CartesianMesh>(
        comm, set_id, block_id, num_i_blocks, num_j_blocks, num_k_blocks,
        global_x_edges.size(), global_y_edges.size(), x_offset, y_offset,
        z_offset, local_x_edges, local_y_edges, local_z_edges );
}

//---------------------------------------------------------------------------//
// Calculate local edge arrays.
void MonteCarloMesh::computeLocalEdges( const std::vector<double> &global_edges,
                                        const std::vector<double> &bnd_mesh,
                                        const int my_block,
                                        std::vector<double> &local_edges,
                                        int &offset ) const
{
    // Find the starting node of the local block.
    auto block_lower = std::lower_bound(
        global_edges.begin(), global_edges.end(), bnd_mesh[my_block] );

    // Because we used lower bound, if the starting node is not the first node
    // then move the starting node back one so we capture the cell in which
    // the boundary lies. If it was the first node, then do not change it as
    // we will start with the first node.
    if ( block_lower != global_edges.begin() )
        --block_lower;

    // The starting node offset is the distance from the front of the edge
    // array to the front of the first cell in the local block.
    offset = std::distance( global_edges.begin(), block_lower );
    DTK_CHECK( offset >= 0 );

    // Find the ending node of the local block. If the boundary mesh is past
    // the end of the global array, this will return the end of that array.
    DTK_CHECK( (std::size_t)my_block + 1 < bnd_mesh.size() );
    auto block_upper = std::upper_bound(
        global_edges.begin(), global_edges.end(), bnd_mesh[my_block + 1] );

    // Because we used upper bound, if the ending node is not the last node
    // then move the ending node back one so we capture the cell in which the
    // boundary lies. If it was the last node then do not change it as we will
    // end with the last node.
    if ( block_upper != global_edges.end() )
        ++block_upper;

    // Calculate the size of the local partition.
    std::size_t local_size = std::distance( block_lower, block_upper );
    DTK_CHECK( local_size <= global_edges.size() );

    // Compose the local edge array.
    local_edges.resize( local_size );
    std::copy( block_lower, block_upper, local_edges.begin() );
}

//---------------------------------------------------------------------------//

} // end namespace Benchmark
} // end namespace DataTransferKit
