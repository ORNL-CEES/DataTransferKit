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
/*!
 * \file DTK_Benchmark_DeterministicMesh.cpp
 * \brief Deterministic mesh interface for the hybrid transport benchmark.
 */
//---------------------------------------------------------------------------//

#include "DTK_DBC.hpp"

#include <algorithm>

namespace DataTransferKit
{
namespace Benchmark
{
//---------------------------------------------------------------------------//
// Uniform cell size constructor. Will determine an "optimal" number of blocks
// based on the communicator size.
DeterministicMesh::DeterministicMesh(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm, const int num_cells_i,
    const int num_cells_j, const int num_cells_k, const double delta_x,
    const double delta_y, const double delta_z )
{
    // Calculate the number of I/J blocks that gives the most "square"
    // decomposition.
    int num_blocks = comm->getSize();
    int num_i_blocks;
    int num_j_blocks;
    num_i_blocks = std::sqrt( static_cast<double>( num_blocks ) );
    num_j_blocks = num_blocks / num_i_blocks;
    while ( num_i_blocks * num_j_blocks != num_blocks )
    {
        DTK_CHECK( num_i_blocks < num_blocks );
        num_i_blocks++;
        num_j_blocks = num_blocks / num_i_blocks;
    }

    // Create uniform global edges and parititon the mesh.
    createEdgesAndPartition( comm, num_i_blocks, num_j_blocks, num_cells_i,
                             num_cells_j, num_cells_k, delta_x, delta_y,
                             delta_z );
}

//---------------------------------------------------------------------------//
// Uniform cell size/block size constructor.
DeterministicMesh::DeterministicMesh(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm, const int num_i_blocks,
    const int num_j_blocks, const int num_cells_i, const int num_cells_j,
    const int num_cells_k, const double delta_x, const double delta_y,
    const double delta_z )
{
    // Create uniform global edges and parititon the mesh.
    createEdgesAndPartition( comm, num_i_blocks, num_j_blocks, num_cells_i,
                             num_cells_j, num_cells_k, delta_x, delta_y,
                             delta_z );
}

//---------------------------------------------------------------------------//
// Global edge/block size constructor.
DeterministicMesh::DeterministicMesh(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm, const int num_i_blocks,
    const int num_j_blocks, const std::vector<double> &global_x_edges,
    const std::vector<double> &global_y_edges,
    const std::vector<double> &global_z_edges )
{
    // Partition the mesh.
    partition( comm, num_i_blocks, num_j_blocks, global_x_edges, global_y_edges,
               global_z_edges );
}

//---------------------------------------------------------------------------//
// Create edges and partition the mesh.
void DeterministicMesh::createEdgesAndPartition(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm, const int num_i_blocks,
    const int num_j_blocks, const int num_cells_i, const int num_cells_j,
    const int num_cells_k, const double delta_x, const double delta_y,
    const double delta_z )
{
    // Create uniform global edge arrays.
    std::vector<double> global_x_edges( num_cells_i + 1 );
    for ( int n = 0; n < num_cells_i + 1; ++n )
        global_x_edges = n * delta_x;
    std::vector<double> global_y_edges( num_cells_j + 1 );
    for ( int n = 0; n < num_cells_j + 1; ++n )
        global_y_edges = n * delta_y;
    std::vector<double> global_z_edges( num_cells_k + 1 );
    for ( int n = 0; n < num_cells_k + 1; ++n )
        global_z_edges = n * delta_z;

    // Partition the mesh.
    partition( comm, num_i_blocks, num_j_blocks, global_x_edges, global_y_edges,
               global_z_edges );
}

//---------------------------------------------------------------------------//
// Partition the mesh.
void DeterministicMesh::partition(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm, const int num_i_blocks,
    const int num_j_blocks, const std::vector<double> &global_x_edges,
    const std::vector<double> &global_y_edges,
    const std::vector<double> &global_z_edges )
{
    // There is only 1 set so set id is 0.
    int set_id = 0;

    // Block id will be the rank in the communicator.
    int block_id = comm->getRank();

    // Compute total number of blocks.
    int num_block = num_i_blocks * num_j_blocks;

    // Determine the ij indices of this block.
    int j_block = block_id / num_i_blocks;
    int i_block = block_id - j_block * num_i_blocks;

    // Get the global cell sizes.
    int x_global_num_cell = global_x_edges - 1;
    int y_global_num_cell = global_y_edges - 1;
    int z_global_num_cell = global_z_edges - 1;

    // This function computes the local number of cells in a given direction
    // in a given block. X and Y will be partitioned, Z will not. Cells are
    // evenly divided amongst blocks. If there are an uneven number of cells
    // per block then each block gets an extra cell until all the cells have
    // been divided amongst blocks.
    auto compute_local_num_cell = []( const int global_num_cell,
                                      const int num_block, const int block ) {
        int local_num = global_num_cell / num_block;
        int append_cells = global_num_cell % num_block;
        if ( block < append_cells )
            ++local_num;
        return local_num;
    };

    // Calculate the offsets into the global edge arrays for all blocks. Do
    // this by looping through all blocks until we get to the current block
    // and adding the number of cells they will have.
    auto compute_offset = []( const int global_num_cell, const int num_block,
                              const int block ) {
        int offset = 0;
        for ( int n = 0; n < block; ++n )
            offset += compute_local_num_cell( global_num_cell, num_block, n );
        return offset;
    };
    int i_offset = compute_offset( x_global_num_cell, num_i_blocks, i_block );
    int j_offset = compute_offset( y_global_num_cell, num_j_blocks, j_block );
    int k_offset = 0;

    // Calculate the local number of cells in each direction for this block.
    int x_local_num_cell =
        compute_local_num_cell( x_global_num_cell, num_i_blocks, i_block );
    int y_local_num_cell =
        compute_local_num_cell( y_global_num_cell, num_j_blocks, j_block );
    int z_local_num_cell = z_global_num_cell;

    // Create the local edge arrays in X and Y.
    std::vector<double> local_x_edges( x_local_num_cell );
    std::copy( global_x_edges.begin() + i_offset,
               global_x_edges.begin() + i_offset + x_local_num_cell,
               local_x_edges.begin() );
    std::vector<double> local_y_edges( y_local_num_cell );
    std::copy( global_y_edges.begin() + j_offset,
               global_y_edges.begin() + j_offset + y_local_num_cell,
               local_y_edges.begin() );

    // Build the mesh data. Note that Z is not partitioned so there is no
    // offset. The global z edge array is not partitioned in this case so just
    // reuse that.
    this->buildMeshData( comm, set_id, block_id, global_x_edges.size(),
                         global_y_edges.size(), global_z_edges.size(), i_offset,
                         j_offset, k_offset, local_x_edges, local_y_edges,
                         global_z_edges );
}

//---------------------------------------------------------------------------//

} // end namespace Benchmark
} // end namespace DataTransferKit
