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

#include "DTK_Benchmark_CartesianMesh.hpp"

#include <Kokkos_Core.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <exception>
#include <vector>

//---------------------------------------------------------------------------//
// Create a cartesian mesh and check the results.
TEUCHOS_UNIT_TEST( CartesianMesh, cartesian_mesh )
{
    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Make sure that the test criteria is satisfied for the comm.
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();
    bool comm_is_good = ( 1 == comm_size || 2 == comm_size || 4 == comm_size ||
                          8 == comm_size );
    if ( !comm_is_good )
        throw std::runtime_error( "Wrong communicator size. Only communicators "
                                  "of size 1, 2, 4, and 8 are valid for this "
                                  "test." );

    // Set the decomposition parameters.
    int set_id = 0;
    int block_id = comm_rank;
    int num_i_blocks = -1;
    int num_j_blocks = -1;
    int num_k_blocks = -1;

    // Set the local number of nodes.
    int x_local_num_node = 12;
    int y_local_num_node = 14;
    int z_local_num_node = 13;
    int x_local_num_cell = x_local_num_node - 1;
    int y_local_num_cell = y_local_num_node - 1;
    int z_local_num_cell = z_local_num_node - 1;
    int local_num_node = x_local_num_node * y_local_num_node * z_local_num_node;

    // Set the global "partitioning" parameters. The transport partitioners
    // can really only be partitioned onto communicators with an even number
    // of processors. This test is written for 1, 2, 4, and 8 processors.
    //
    // The offsets here indicate the starting index of this rank's nodes in
    // the "global" node arrays. The cartesian mesh only needs the local node
    // arrays but we need these offsets to compose global indices.
    int x_global_num_cell = 0;
    int y_global_num_cell = 0;
    int x_offset = 0;
    int y_offset = 0;
    int z_offset = 0;

    // Comm size 1 has no partitioning
    /*
      +--------------+
      |              |
      |              |
      |      0       |
      |              |
      |              |
      +--------------+
     */
    if ( 1 == comm_size )
    {
        x_global_num_cell = x_local_num_cell;
        y_global_num_cell = y_local_num_cell;

        x_offset = 0;
        y_offset = 0;
        z_offset = 0;

        num_i_blocks = 1;
        num_j_blocks = 1;
        num_k_blocks = 1;
    }

    // Comm size 2 has partitioning in x
    /*
      +--------------+--------------+
      |              |              |
      |              |              |    y
      |      0       |      1       |    ^
      |              |              |    |
      |              |              |    |
      +--------------+--------------+    ---->x
     */
    else if ( 2 == comm_size )
    {
        x_global_num_cell = 2 * x_local_num_cell;
        y_global_num_cell = y_local_num_cell;

        // Rank 1 gets an x offset in this case.
        x_offset = ( 1 == comm_rank ) ? x_local_num_cell : 0;

        // No offsets in y an z
        y_offset = 0;
        z_offset = 0;

        num_i_blocks = 2;
        num_j_blocks = 1;
    }

    // Comm size 4 has partitioning in x and y
    /*
      +--------------+--------------+
      |              |              |
      |              |              |
      |      2       |      3       |
      |              |              |
      |              |              |
      +--------------+--------------+
      |              |              |
      |              |              |    y
      |      0       |      1       |    ^
      |              |              |    |
      |              |              |    |
      +--------------+--------------+    ---->x
     */
    else if ( 4 == comm_size )
    {
        x_global_num_cell = 2 * x_local_num_cell;
        y_global_num_cell = 2 * y_local_num_cell;

        // Ranks get 1 and 3 get an x offset
        x_offset = ( 1 == comm_rank || 3 == comm_rank ) ? x_local_num_cell : 0;

        // Ranks 2 and 3 get a y offset.
        y_offset = ( 2 == comm_rank || 3 == comm_rank ) ? y_local_num_cell : 0;

        // no z offset
        z_offset = 0;

        num_i_blocks = 2;
        num_j_blocks = 2;
        num_k_blocks = 1;
    }

    // Comm size 8 has partitioning in x, y, and z (rank 2 is hidden in the
    // back lower left corner of the drawing).
    /*
           +------------+-----------+
          /            /           /|
         /     6      /     7     / |
        /            /           /  |
        +------------+-----------+  |
       /            /           /|7 +
      /     4      /    5      / |  /
     /            /           /  | /|
     +------------+-----------+  |/ |
     |            |           |5 +  |
     |            |           |  / 3|
     |     4      |     5     | / | +
     |            |           |/  | /
     +------------+-----------+   |/     z    y
     |            |           | 1 +      ^   ~
     |            |           |  /       |  /
     |     0      |     1     | /        | /
     |            |           |/         |/
     +------------+-----------+          ----->x
     */
    else if ( 8 == comm_size )
    {
        x_global_num_cell = 2 * x_local_num_cell;
        y_global_num_cell = 2 * y_local_num_cell;

        // Ranks 1, 3, 5 and 7 get an x offset.
        x_offset = ( 1 == comm_rank || 3 == comm_rank || 5 == comm_rank ||
                     7 == comm_rank )
                       ? x_local_num_cell
                       : 0;

        // Ranks 2, 3, 6, and 7 get a y offset.
        y_offset = ( 2 == comm_rank || 3 == comm_rank || 6 == comm_rank ||
                     7 == comm_rank )
                       ? y_local_num_cell
                       : 0;

        // Ranks 4, 5, 6, and 7 get a z offset.
        z_offset = ( 4 == comm_rank || 5 == comm_rank || 6 == comm_rank ||
                     7 == comm_rank )
                       ? z_local_num_cell
                       : 0;

        num_i_blocks = 2;
        num_j_blocks = 2;
        num_k_blocks = 2;
    }

    // Set the local number of cells.
    int local_num_cell = x_local_num_cell * y_local_num_cell * z_local_num_cell;

    // Set the global number of cells.
    int x_global_num_node = x_global_num_cell + 1;
    int y_global_num_node = y_global_num_cell + 1;

    // Create local edges. The global low corner of the mesh will be
    // (0,0,0). Use a mesh spacing parameter in each direction to build the
    // edges.
    std::vector<double> local_x_edges( x_local_num_node );
    std::vector<double> local_y_edges( y_local_num_node );
    std::vector<double> local_z_edges( z_local_num_node );
    double dx = 0.5;
    double dy = 0.4;
    double dz = 0.3;
    double x_start = dx * x_offset;
    double y_start = dy * y_offset;
    double z_start = dz * z_offset;
    for ( int n = 0; n < x_local_num_node; ++n )
        local_x_edges[n] = x_start + n * dx;
    for ( int n = 0; n < y_local_num_node; ++n )
        local_y_edges[n] = y_start + n * dy;
    for ( int n = 0; n < z_local_num_node; ++n )
        local_z_edges[n] = z_start + n * dz;

    // Create the mesh.
    DataTransferKit::Benchmark::CartesianMesh mesh(
        comm, set_id, block_id, num_i_blocks, num_j_blocks, num_k_blocks,
        x_global_num_node, y_global_num_node, x_offset, y_offset, z_offset,
        local_x_edges, local_y_edges, local_z_edges );

    // Check the mesh decomposition parameters.
    TEST_EQUALITY( mesh.setId(), set_id );
    TEST_EQUALITY( mesh.blockId(), block_id );

    // Check the global ids of the nodes with an ijk indexer.
    auto node_ids = mesh.localNodeGlobalIds();
    auto node_ids_host = Kokkos::create_mirror_view( node_ids );
    Kokkos::deep_copy( node_ids_host, node_ids );
    TEST_EQUALITY( Teuchos::as<int>( node_ids_host.extent( 0 ) ),
                   local_num_node );
    auto local_node_id = [=]( const int i, const int j, const int k ) {
        return i + j * x_local_num_node +
               k * x_local_num_node * y_local_num_node;
    };
    auto global_node_id = [=]( const int i, const int j, const int k ) {
        return ( i + x_offset ) + ( j + y_offset ) * x_global_num_node +
               ( k + z_offset ) * x_global_num_node * y_global_num_node;
    };
    for ( int k = 0; k < z_local_num_node; ++k )
    {
        for ( int j = 0; j < y_local_num_node; ++j )
        {
            for ( int i = 0; i < x_local_num_node; ++i )
            {
                auto lid = local_node_id( i, j, k );
                auto gid = global_node_id( i, j, k );
                TEST_EQUALITY( Teuchos::as<int>( node_ids_host( lid ) ), gid );
            }
        }
    }

    // Check the local coordinates of the nodes.
    auto node_coords = mesh.localNodeCoordinates();
    auto node_coords_host = Kokkos::create_mirror_view( node_coords );
    Kokkos::deep_copy( node_coords_host, node_coords );
    for ( int k = 0; k < z_local_num_node; ++k )
    {
        for ( int j = 0; j < y_local_num_node; ++j )
        {
            for ( int i = 0; i < x_local_num_node; ++i )
            {
                auto node_id = local_node_id( i, j, k );

                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( node_id, 0 ) ),
                    x_start + i * dx, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( node_id, 1 ) ),
                    y_start + j * dy, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( node_id, 2 ) ),
                    z_start + k * dz, 1e-14 );
            }
        }
    }

    // Check the global ids of the cells.
    auto cell_ids = mesh.localCellGlobalIds();
    auto cell_ids_host = Kokkos::create_mirror_view( cell_ids );
    Kokkos::deep_copy( cell_ids_host, cell_ids );
    TEST_EQUALITY( Teuchos::as<int>( cell_ids.extent( 0 ) ), local_num_cell );
    auto local_cell_id = [=]( const int i, const int j, const int k ) {
        return i + j * x_local_num_cell +
               k * x_local_num_cell * y_local_num_cell;
    };
    auto global_cell_id = [=]( const int i, const int j, const int k ) {
        return ( i + x_offset ) + ( j + y_offset ) * x_global_num_cell +
               ( k + z_offset ) * x_global_num_cell * y_global_num_cell;
    };
    for ( int k = 0; k < z_local_num_cell; ++k )
    {
        for ( int j = 0; j < y_local_num_cell; ++j )
        {
            for ( int i = 0; i < x_local_num_cell; ++i )
            {
                auto lid = local_cell_id( i, j, k );
                auto gid = global_cell_id( i, j, k );
                TEST_EQUALITY( Teuchos::as<int>( cell_ids_host( lid ) ), gid );
            }
        }
    }

    // Check the coordinates of the cell centers and the connectivity of the
    // cells.
    auto cell_conn = mesh.localCellConnectivity();
    auto cell_conn_host = Kokkos::create_mirror_view( cell_conn );
    Kokkos::deep_copy( cell_conn_host, cell_conn );
    auto cell_coords = mesh.localCellCenterCoordinates();
    auto cell_coords_host = Kokkos::create_mirror_view( cell_coords );
    Kokkos::deep_copy( cell_coords_host, cell_coords );
    for ( int k = 0; k < z_local_num_cell; ++k )
    {
        for ( int j = 0; j < y_local_num_cell; ++j )
        {
            for ( int i = 0; i < x_local_num_cell; ++i )
            {
                auto cell_id = local_cell_id( i, j, k );

                // Check the center of the cell.
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( cell_coords_host( cell_id, 0 ) ),
                    x_start + ( i + 0.5 ) * dx, 1.0e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( cell_coords_host( cell_id, 1 ) ),
                    y_start + ( j + 0.5 ) * dy, 1.0e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( cell_coords_host( cell_id, 2 ) ),
                    z_start + ( k + 0.5 ) * dz, 1.0e-14 );

                // Check the connectivity by checking the coordinates. This
                // lets us use the connectivity array to index back into the
                // node coordinates to make sure we really constructed the
                // correct cells.
                auto cell_node = cell_conn_host( cell_id, 0 );
                int ni = i;
                int nj = j;
                int nk = k;
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 0 ) ),
                    x_start + ni * dx, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 1 ) ),
                    y_start + nj * dy, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 2 ) ),
                    z_start + nk * dz, 1e-14 );

                cell_node = cell_conn_host( cell_id, 1 );
                ni = i + 1;
                nj = j;
                nk = k;
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 0 ) ),
                    x_start + ni * dx, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 1 ) ),
                    y_start + nj * dy, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 2 ) ),
                    z_start + nk * dz, 1e-14 );

                cell_node = cell_conn_host( cell_id, 2 );
                ni = i + 1;
                nj = j + 1;
                nk = k;
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 0 ) ),
                    x_start + ni * dx, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 1 ) ),
                    y_start + nj * dy, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 2 ) ),
                    z_start + nk * dz, 1e-14 );

                cell_node = cell_conn_host( cell_id, 3 );
                ni = i;
                nj = j + 1;
                nk = k;
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 0 ) ),
                    x_start + ni * dx, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 1 ) ),
                    y_start + nj * dy, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 2 ) ),
                    z_start + nk * dz, 1e-14 );

                cell_node = cell_conn_host( cell_id, 4 );
                ni = i;
                nj = j;
                nk = k + 1;
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 0 ) ),
                    x_start + ni * dx, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 1 ) ),
                    y_start + nj * dy, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 2 ) ),
                    z_start + nk * dz, 1e-14 );

                cell_node = cell_conn_host( cell_id, 5 );
                ni = i + 1;
                nj = j;
                nk = k + 1;
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 0 ) ),
                    x_start + ni * dx, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 1 ) ),
                    y_start + nj * dy, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 2 ) ),
                    z_start + nk * dz, 1e-14 );

                cell_node = cell_conn_host( cell_id, 6 );
                ni = i + 1;
                nj = j + 1;
                nk = k + 1;
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 0 ) ),
                    x_start + ni * dx, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 1 ) ),
                    y_start + nj * dy, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 2 ) ),
                    z_start + nk * dz, 1e-14 );

                cell_node = cell_conn_host( cell_id, 7 );
                ni = i;
                nj = j + 1;
                nk = k + 1;
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 0 ) ),
                    x_start + ni * dx, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 1 ) ),
                    y_start + nj * dy, 1e-14 );
                TEST_FLOATING_EQUALITY(
                    static_cast<double>( node_coords_host( cell_node, 2 ) ),
                    z_start + nk * dz, 1e-14 );
            }
        }
    }
}

//---------------------------------------------------------------------------//
