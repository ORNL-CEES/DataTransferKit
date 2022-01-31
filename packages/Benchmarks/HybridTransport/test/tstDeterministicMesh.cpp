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

#include "DTK_Benchmark_DeterministicMesh.hpp"

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

    // Cell counts;
    int x_global_num_cell;
    int y_global_num_cell;
    int z_global_num_cell;

    // Cell sizes.
    double dx;
    double dy;
    double dz;
};

//---------------------------------------------------------------------------//
// Create the parallel problem.
TestParameters createProblem()
{
    // Parameters.
    TestParameters p;

    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Make sure that the test criteria is satisfied for the comm.
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();
    bool comm_is_good = ( 1 == comm_size || 0 == ( comm_size % 2 ) );
    if ( !comm_is_good )
        throw std::runtime_error(
            "Wrong communicator size. Only communicators "
            "of size 1 or of even size are valid for this "
            "test." );

    // Set the decomposition parameters.
    p.set_id = 0;
    p.block_id = comm_rank;

    // Set the global number of cells.
    p.x_global_num_cell = 7 * comm_size;
    p.y_global_num_cell = 8 * comm_size;
    p.z_global_num_cell = 6 * comm_size;

    // Set the cell size.
    p.dx = 0.4;
    p.dy = 0.5;
    p.dz = 0.3;

    // Return the parameters.
    return p;
}

//---------------------------------------------------------------------------//
// Test the results.
void checkResults( const DataTransferKit::Benchmark::CartesianMesh &mesh,
                   const TestParameters &p, bool &success,
                   Teuchos::FancyOStream &out )
{
    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();

    // Check the mesh decomposition parameters.
    TEST_EQUALITY( mesh.setId(), p.set_id );
    TEST_EQUALITY( mesh.blockId(), p.block_id );

    // Check that we created the right number of sets and blocks.
    TEST_EQUALITY( 1, mesh.numSets() );
    TEST_EQUALITY( comm_size, mesh.numBlocks() );

    // Check that the local cell partitions are of the right size (i.e. the
    // local cell sizes add up to the global mesh sizes).
    auto cell_ids = mesh.localCellGlobalIds();
    int local_num_cell = cell_ids.extent( 0 );
    int global_num_cell = 0;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, local_num_cell,
                        Teuchos::ptrFromRef( global_num_cell ) );
    int expected_num_cell =
        p.x_global_num_cell * p.y_global_num_cell * p.z_global_num_cell;
    TEST_EQUALITY( global_num_cell, expected_num_cell );

    // Also make sure that the local partition has at least 1 cell.
    TEST_ASSERT( 0 < local_num_cell );

    // Check the actual distribution of cells by computing the sum of the
    // global cell ids owned by each rank. This sum should be unique and
    // correct relative to a reference only if all cells were distributed
    // properly to unique owning processors. In other words, each cell id
    // should contribute to the global sum once and only once.
    //
    // This also checks the offset computation of the partitioner because the
    // offsets are used to calculate global ids. If the offsets are wrong, the
    // global ids are wrong, and therefore the sum will be wrong. We know that
    // the id computation is correct because we already tested that when we
    // tested CartesianMesh
    std::size_t expected_id_sum = 0;
    for ( int i = 0; i < expected_num_cell; ++i )
        expected_id_sum += i;
    std::size_t local_id_sum = 0;
    auto cell_ids_host = Kokkos::create_mirror_view( cell_ids );
    Kokkos::deep_copy( cell_ids_host, cell_ids );
    for ( int i = 0; i < local_num_cell; ++i )
        local_id_sum += cell_ids_host( i );
    std::size_t global_id_sum = 0;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, local_id_sum,
                        Teuchos::ptrFromRef( global_id_sum ) );
    TEST_EQUALITY( global_id_sum, expected_id_sum );
}

//---------------------------------------------------------------------------//
// Create a deterministic mesh with uniform cells.
TEUCHOS_UNIT_TEST( DeterministicMesh, uniform_cell_constructor )
{
    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Create test parameters.
    auto params = createProblem();

    // Build a deterministic mesh.
    DataTransferKit::Benchmark::DeterministicMesh mesh(
        comm, params.x_global_num_cell, params.y_global_num_cell,
        params.z_global_num_cell, params.dx, params.dy, params.dz );

    // Check the result.
    checkResults( *( mesh.cartesianMesh() ), params, success, out );
}

//---------------------------------------------------------------------------//
// Create a deterministic mesh with the global edge constructor.
TEUCHOS_UNIT_TEST( DeterministicMesh, global_edge_constructor )
{
    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Create test parameters.
    auto params = createProblem();

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
    DataTransferKit::Benchmark::DeterministicMesh mesh(
        comm, global_x_edges, global_y_edges, global_z_edges );

    // Check the result.
    checkResults( *( mesh.cartesianMesh() ), params, success, out );
}

//---------------------------------------------------------------------------//
