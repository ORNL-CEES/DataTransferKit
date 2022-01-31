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

#include "DTK_Benchmark_MonteCarloMesh.hpp"

#include <Kokkos_Core.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <exception>
#include <vector>

//---------------------------------------------------------------------------//
// Test paramters.
struct TestParameters
{
    // Parallel decomposition.
    int set_id;
    int block_id;
    int num_sets;
    int num_blocks;

    // Cell counts;
    int x_global_num_cell;
    int y_global_num_cell;
    int z_global_num_cell;

    // Cell sizes.
    double dx;
    double dy;
    double dz;

    // Boundary mesh.
    std::vector<double> x_bnd_mesh;
    std::vector<double> y_bnd_mesh;
    std::vector<double> z_bnd_mesh;
};

//---------------------------------------------------------------------------//
// Create the parallel problem.
TestParameters createProblem( const int num_sets )
{
    // Parameters.
    TestParameters p;

    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Make sure that the test criteria is satisfied for the comm.
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();
    bool comm_is_good = ( 1 == comm_size || 2 == comm_size || 4 == comm_size ||
                          8 == comm_size );
    if ( !comm_is_good )
        throw std::runtime_error(
            "Wrong communicator size. Only communicators "
            "of size 1, 2, 4, or 8 are valid for this test." );

    // Number of blocks and sets
    p.num_sets = num_sets;
    p.num_blocks = comm_size / num_sets;

    // Set the decomposition parameters.
    p.set_id = std::floor( comm_rank / p.num_blocks );
    p.block_id = comm_rank - ( p.set_id * p.num_blocks );

    // Set the global number of cells.
    p.x_global_num_cell = 24;
    p.y_global_num_cell = 16;
    p.z_global_num_cell = 20;

    // Set the cell size.
    p.dx = 0.4;
    p.dy = 0.5;
    p.dz = 0.3;

    // Create the boundary mesh.
    double x_max = p.x_global_num_cell * p.dx;
    double y_max = p.y_global_num_cell * p.dy;
    double z_max = p.z_global_num_cell * p.dz;
    if ( 1 == p.num_blocks )
    {
        p.x_bnd_mesh = {-0.1, x_max + 0.1};
        p.y_bnd_mesh = {-0.1, y_max + 0.1};
        p.z_bnd_mesh = {-0.1, z_max + 0.1};
    }
    else if ( 2 == p.num_blocks )
    {
        p.x_bnd_mesh = {-0.1, x_max / 1.93, x_max + 0.1};
        p.y_bnd_mesh = {-0.1, y_max + 0.1};
        p.z_bnd_mesh = {-0.1, z_max + 0.1};
    }
    else if ( 4 == p.num_blocks )
    {
        p.x_bnd_mesh = {-0.1, x_max / 1.93, x_max + 0.1};
        p.y_bnd_mesh = {-0.1, y_max / 1.93, y_max + 0.1};
        p.z_bnd_mesh = {-0.1, z_max + 0.1};
    }
    else if ( 8 == p.num_blocks )
    {
        p.x_bnd_mesh = {-0.1, x_max / 1.93, x_max + 0.1};
        p.y_bnd_mesh = {-0.1, y_max / 1.93, y_max + 0.1};
        p.z_bnd_mesh = {-0.1, z_max / 1.93, z_max + 0.1};
    }

    // Return the parameters.
    return p;
}

//---------------------------------------------------------------------------//
// Test the results.
void checkResults( const int num_sets,
                   const DataTransferKit::Benchmark::CartesianMesh &mesh,
                   const TestParameters &p, bool &success,
                   Teuchos::FancyOStream &out )
{
    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();

    // Check the mesh decomposition parameters.
    TEST_EQUALITY( mesh.setId(), p.set_id );
    TEST_EQUALITY( mesh.blockId(), p.block_id );

    // Check that we created the right number of blocks and sets.
    TEST_EQUALITY( num_sets, mesh.numSets() );
    TEST_EQUALITY( comm_size, mesh.numSets() * mesh.numBlocks() );

    // Check that the local cell partitions are of the right size (i.e. the
    // local cell sizes add up to the global mesh sizes plus some overlap due
    // to ghosting from the intersection with the boundary mesh).

    // Start by getting the local number of cells in the mesh and then summing
    // to get the global unique+ghosted number.
    auto cell_ids = mesh.localCellGlobalIds();
    int local_num_cell = cell_ids.extent( 0 );
    int global_num_cell = 0;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, local_num_cell,
                        Teuchos::ptrFromRef( global_num_cell ) );

    // Now calculate the expected number of cells. First get the uniquely
    // owned number.
    int expected_num_cell = num_sets * p.x_global_num_cell *
                            p.y_global_num_cell * p.z_global_num_cell;

    // Then add ghosting.
    expected_num_cell += ( p.num_blocks > 1 ) ? num_sets * p.y_global_num_cell *
                                                    p.z_global_num_cell
                                              : 0;
    expected_num_cell +=
        ( p.num_blocks > 2 )
            ? num_sets * ( p.x_global_num_cell + 1 ) * p.z_global_num_cell
            : 0;
    expected_num_cell += ( p.num_blocks > 4 )
                             ? num_sets * ( p.x_global_num_cell + 1 ) *
                                   ( p.y_global_num_cell + 1 )
                             : 0;
    TEST_EQUALITY( global_num_cell, expected_num_cell );

    // Also make sure that the local partition has at least 1 cell.
    TEST_ASSERT( 0 < local_num_cell );

    // Check the actual distribution of cells by computing the sum of the
    // global cell ids owned by each rank. This sum should be unique and
    // correct relative to a reference only if all cells were distributed
    // properly to unique and ghosted owning processors. In other words, each
    // cell id should contribute to the global sum once per set if it is not
    // on a partition boundary and twice per set if it is.
    //
    // This also checks the offset computation of the partitioner because the
    // offsets are used to calculate global ids. If the offsets are wrong, the
    // global ids are wrong, and therefore the sum will be wrong. We know that
    // the id computation is correct because we already tested that when we
    // tested CartesianMesh

    // Then compute the expected sum without ghosting.
    std::size_t expected_id_sum = 0;
    int global_num_unique_cell =
        p.x_global_num_cell * p.y_global_num_cell * p.z_global_num_cell;
    for ( int i = 0; i < global_num_unique_cell; ++i )
        expected_id_sum += num_sets * i;

    // Then add the ids (if any) that are double counted.
    if ( p.num_blocks > 1 )
    {
        // There is one layer in X that is double counted due to overlap.
        int x_id = std::floor( p.x_bnd_mesh[1] / p.dx );
        for ( int j = 0; j < p.y_global_num_cell; ++j )
            for ( int k = 0; k < p.z_global_num_cell; ++k )
                expected_id_sum += num_sets * ( x_id + p.x_global_num_cell * j +
                                                p.x_global_num_cell *
                                                    p.y_global_num_cell * k );
    }
    if ( p.num_blocks > 2 )
    {
        // There is another layer in Y that is additionaly double counted.
        int y_id = std::floor( p.y_bnd_mesh[1] / p.dy );
        for ( int i = 0; i < p.x_global_num_cell; ++i )
            for ( int k = 0; k < p.z_global_num_cell; ++k )
                expected_id_sum += num_sets * ( i + p.x_global_num_cell * y_id +
                                                p.x_global_num_cell *
                                                    p.y_global_num_cell * k );

        // And a a second row of cells in Z that is in both overlapping zones
        // and therefore counted yet another time.
        int x_id = std::floor( p.x_bnd_mesh[1] / p.dx );
        for ( int k = 0; k < p.z_global_num_cell; ++k )
            expected_id_sum +=
                num_sets * ( x_id + p.x_global_num_cell * y_id +
                             p.x_global_num_cell * p.y_global_num_cell * k );
    }
    if ( p.num_blocks > 4 )
    {
        // There is another layer in Z that is additionally double counted.
        int z_id = std::floor( p.z_bnd_mesh[1] / p.dz );
        for ( int i = 0; i < p.x_global_num_cell; ++i )
            for ( int j = 0; j < p.y_global_num_cell; ++j )
                expected_id_sum +=
                    num_sets *
                    ( i + p.x_global_num_cell * j +
                      p.x_global_num_cell * p.y_global_num_cell * z_id );

        // And a a second row of cells in X that is in both overlapping zones
        // and therefore counted yet another time.
        int x_id = std::floor( p.x_bnd_mesh[1] / p.dx );
        for ( int j = 0; j < p.y_global_num_cell; ++j )
            expected_id_sum +=
                num_sets * ( x_id + p.x_global_num_cell * j +
                             p.x_global_num_cell * p.y_global_num_cell * z_id );

        // And a a second row of cells in Y that is in both overlapping zones
        // and therefore counted yet another time.
        int y_id = std::floor( p.y_bnd_mesh[1] / p.dy );
        for ( int i = 0; i < p.x_global_num_cell; ++i )
            expected_id_sum +=
                num_sets * ( i + p.x_global_num_cell * y_id +
                             p.x_global_num_cell * p.y_global_num_cell * z_id );

        // And one cell that intersects all overlapping regions between blocks.
        expected_id_sum +=
            num_sets * ( x_id + p.x_global_num_cell * y_id +
                         p.x_global_num_cell * p.y_global_num_cell * z_id );
    }

    // Then get the actual sum.
    std::size_t local_id_sum = 0;
    auto cell_ids_host = Kokkos::create_mirror_view( cell_ids );
    Kokkos::deep_copy( cell_ids_host, cell_ids );
    for ( int i = 0; i < local_num_cell; ++i )
        local_id_sum += cell_ids_host( i );
    std::size_t global_id_sum = 0;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, local_id_sum,
                        Teuchos::ptrFromRef( global_id_sum ) );

    // Compare.
    TEST_EQUALITY( global_id_sum, expected_id_sum );
}

//---------------------------------------------------------------------------//
// Run a uniform cell test.
void uniformCellTest( const int num_sets, bool &success,
                      Teuchos::FancyOStream &out )
{
    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Only run this test if we have enough cores for the sets.
    if ( num_sets <= comm->getSize() )
    {
        // Create test parameters.
        auto params = createProblem( num_sets );

        // Build a deterministic mesh.
        DataTransferKit::Benchmark::MonteCarloMesh mesh(
            comm, params.num_sets, params.x_global_num_cell,
            params.y_global_num_cell, params.z_global_num_cell, params.dx,
            params.dy, params.dz, params.x_bnd_mesh, params.y_bnd_mesh,
            params.z_bnd_mesh );

        // Check the result.
        checkResults( num_sets, *( mesh.cartesianMesh() ), params, success,
                      out );
    }

    else
    {
        std::cout << "Not enough cores to test with so " << num_sets
                  << " sets test not run." << std::endl;
    }
}

//---------------------------------------------------------------------------//
// Create a deterministic mesh with uniform cells and 1 set.
TEUCHOS_UNIT_TEST( MonteCarloMesh, 1_set )
{
    uniformCellTest( 1, success, out );
}

//---------------------------------------------------------------------------//
// Create a deterministic mesh with uniform cells and 2 sets.
TEUCHOS_UNIT_TEST( MonteCarloMesh, 2_set )
{
    uniformCellTest( 2, success, out );
}

//---------------------------------------------------------------------------//
// Create a deterministic mesh with uniform cells and 4 sets.
TEUCHOS_UNIT_TEST( MonteCarloMesh, 4_set )
{
    uniformCellTest( 4, success, out );
}

//---------------------------------------------------------------------------//
// Create a deterministic mesh with the global edge constructor.
TEUCHOS_UNIT_TEST( MonteCarloMesh, global_edge_constructor )
{
    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Create test parameters.
    int num_sets = 1;
    auto params = createProblem( num_sets );

    // Create global edges.
    std::vector<double> global_x_edges( params.x_global_num_cell + 1 );
    for ( int n = 0; n < params.x_global_num_cell + 1; ++n )
        global_x_edges[n] = params.dx * n;
    std::vector<double> global_y_edges( params.y_global_num_cell + 1 );
    for ( int n = 0; n < params.y_global_num_cell + 1; ++n )
        global_y_edges[n] = params.dy * n;
    std::vector<double> global_z_edges( params.z_global_num_cell + 1 );
    for ( int n = 0; n < params.z_global_num_cell + 1; ++n )
        global_z_edges[n] = params.dz * n;

    // Build a deterministic mesh.
    DataTransferKit::Benchmark::MonteCarloMesh mesh(
        comm, params.num_sets, global_x_edges, global_y_edges, global_z_edges,
        params.x_bnd_mesh, params.y_bnd_mesh, params.z_bnd_mesh );

    // Check the result.
    checkResults( num_sets, *( mesh.cartesianMesh() ), params, success, out );
}

//---------------------------------------------------------------------------//
