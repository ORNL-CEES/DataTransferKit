/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
#ifndef DTK_TESTAPPLICATIONHELPERS_HPP
#define DTK_TESTAPPLICATIONHELPERS_HPP

#include <DTK_CellTypes.h>

#include <Teuchos_UnitTestHarness.hpp>

#include "DTK_APIConstants.h"

template <class UserApplication>
void test_node_list( UserApplication &user_app, Teuchos::FancyOStream &out,
                     bool &success )
{
    // Get a node list.
    auto node_list = user_app.getNodeList();

    // Check the node list.
    auto host_coordinates = Kokkos::create_mirror_view( node_list.coordinates );
    Kokkos::deep_copy( host_coordinates, node_list.coordinates );
    for ( unsigned i = 0; i < SIZE_1; ++i )
    {
        for ( unsigned d = 0; d < SPACE_DIM; ++d )
            TEST_EQUALITY( host_coordinates( i, d ), i + d + OFFSET );
    }
}

template <class UserApplication>
void test_bounding_volume_list( UserApplication &user_app,
                                Teuchos::FancyOStream &out, bool &success )
{
    // Get a bounding volume list.
    auto bv_list = user_app.getBoundingVolumeList();

    // Check the bounding volumes.
    auto host_bounding_volumes =
        Kokkos::create_mirror_view( bv_list.bounding_volumes );
    Kokkos::deep_copy( host_bounding_volumes, bv_list.bounding_volumes );
    for ( unsigned i = 0; i < SIZE_1; ++i )
    {
        for ( unsigned d = 0; d < SPACE_DIM; ++d )
            for ( unsigned b = 0; b < 2; ++b )
                TEST_EQUALITY( host_bounding_volumes( i, d, b ),
                               i + d + b + OFFSET );
    }
}

template <class UserApplication>
void test_polyhedron_list( UserApplication &user_app,
                           Teuchos::FancyOStream &out, bool &success )
{
    // Get a polyhedron list.
    auto poly_list = user_app.getPolyhedronList();

    // Check the list.
    auto host_coordinates = Kokkos::create_mirror_view( poly_list.coordinates );
    Kokkos::deep_copy( host_coordinates, poly_list.coordinates );
    auto host_faces = Kokkos::create_mirror_view( poly_list.faces );
    Kokkos::deep_copy( host_faces, poly_list.faces );
    auto host_nodes_per_face =
        Kokkos::create_mirror_view( poly_list.nodes_per_face );
    Kokkos::deep_copy( host_nodes_per_face, poly_list.nodes_per_face );
    auto host_cells = Kokkos::create_mirror_view( poly_list.cells );
    Kokkos::deep_copy( host_cells, poly_list.cells );
    auto host_faces_per_cell =
        Kokkos::create_mirror_view( poly_list.faces_per_cell );
    Kokkos::deep_copy( host_faces_per_cell, poly_list.faces_per_cell );
    auto host_face_orientation =
        Kokkos::create_mirror_view( poly_list.face_orientation );
    Kokkos::deep_copy( host_face_orientation, poly_list.face_orientation );
    for ( unsigned i = 0; i < SIZE_1; ++i )
    {
        for ( unsigned d = 0; d < SPACE_DIM; ++d )
            TEST_EQUALITY( host_coordinates( i, d ), i + d + OFFSET );
        TEST_EQUALITY( host_faces( i ), i + OFFSET );
        TEST_EQUALITY( host_nodes_per_face( i ), i + OFFSET );
        TEST_EQUALITY( host_cells( i ), i + OFFSET );
        TEST_EQUALITY( host_faces_per_cell( i ), i + OFFSET );
        TEST_EQUALITY( host_face_orientation( i ), 1 );
    }
}

template <class UserApplication>
void test_multiple_topology_cell( UserApplication &user_app,
                                  Teuchos::FancyOStream &out, bool &success )
{
    // Get a cell list.
    auto cell_list = user_app.getCellList();

    // Check the list.
    auto host_coordinates = Kokkos::create_mirror_view( cell_list.coordinates );
    Kokkos::deep_copy( host_coordinates, cell_list.coordinates );
    auto host_cells = Kokkos::create_mirror_view( cell_list.cells );
    Kokkos::deep_copy( host_cells, cell_list.cells );
    auto host_cell_topologies =
        Kokkos::create_mirror_view( cell_list.cell_topologies );
    Kokkos::deep_copy( host_cell_topologies, cell_list.cell_topologies );
    for ( unsigned i = 0; i < SIZE_1; ++i )
    {
        for ( unsigned d = 0; d < SPACE_DIM; ++d )
            TEST_EQUALITY( host_coordinates( i, d ), i + d + OFFSET );
        TEST_EQUALITY( host_cells( i ), i + OFFSET );
        TEST_EQUALITY( host_cell_topologies( i ), DTK_TET_4 );
    }
}

template <class UserApplication>
void test_boundary( UserApplication &user_app, Teuchos::FancyOStream &out,
                    bool &success )
{
    // Test with a cell list.
    {
        // Create a cell list.
        auto cell_list = user_app.getCellList();

        // Get the boundary of the list.
        user_app.getBoundary( cell_list );

        // Check the boundary.
        auto host_boundary_cells =
            Kokkos::create_mirror_view( cell_list.boundary_cells );
        Kokkos::deep_copy( host_boundary_cells, cell_list.boundary_cells );
        auto host_cell_faces_on_boundary =
            Kokkos::create_mirror_view( cell_list.cell_faces_on_boundary );
        Kokkos::deep_copy( host_cell_faces_on_boundary,
                           cell_list.cell_faces_on_boundary );
        for ( unsigned i = 0; i < SIZE_1; ++i )
        {
            TEST_EQUALITY( host_boundary_cells( i ), i + OFFSET );
            TEST_EQUALITY( host_cell_faces_on_boundary( i ), i + OFFSET );
        }
    }

    // Test with a polyhedron list.
    {
        // Create a polyhedron list.
        auto poly_list = user_app.getPolyhedronList();

        // Get the boundary of the list.
        user_app.getBoundary( poly_list );

        // Check the boundary.
        auto host_boundary_cells =
            Kokkos::create_mirror_view( poly_list.boundary_cells );
        Kokkos::deep_copy( host_boundary_cells, poly_list.boundary_cells );
        auto host_cell_faces_on_boundary =
            Kokkos::create_mirror_view( poly_list.cell_faces_on_boundary );
        Kokkos::deep_copy( host_cell_faces_on_boundary,
                           poly_list.cell_faces_on_boundary );
        for ( unsigned i = 0; i < SIZE_1; ++i )
        {
            TEST_EQUALITY( host_boundary_cells( i ), i + OFFSET );
            TEST_EQUALITY( host_cell_faces_on_boundary( i ), i + OFFSET );
        }
    }
}

template <class UserApplication>
void test_adjacency_list( UserApplication &user_app, Teuchos::FancyOStream &out,
                          bool &success )
{
    // Test with a cell list.
    {
        // Create a cell list.
        auto cell_list = user_app.getCellList();

        // Get the adjacency of the list.
        user_app.getAdjacencyList( cell_list );

        // Check the adjacency list.
        auto host_cell_global_ids =
            Kokkos::create_mirror_view( cell_list.cell_global_ids );
        Kokkos::deep_copy( host_cell_global_ids, cell_list.cell_global_ids );
        auto host_adjacent_cells =
            Kokkos::create_mirror_view( cell_list.adjacent_cells );
        Kokkos::deep_copy( host_adjacent_cells, cell_list.adjacent_cells );
        auto host_adjacencies_per_cell =
            Kokkos::create_mirror_view( cell_list.adjacencies_per_cell );
        Kokkos::deep_copy( host_adjacencies_per_cell,
                           cell_list.adjacencies_per_cell );
        for ( unsigned i = 0; i < SIZE_1; ++i )
        {
            TEST_EQUALITY( host_cell_global_ids( i ), i + OFFSET );
            TEST_EQUALITY( host_adjacent_cells( i ), i );
            TEST_EQUALITY( host_adjacencies_per_cell( i ), 1 );
        }
    }

    // Test with a polyhedron list.
    {
        // Create a polyhedron list.
        auto poly_list = user_app.getPolyhedronList();

        // Get the adjacency of the list.
        user_app.getAdjacencyList( poly_list );

        // Check the adjacency list.
        auto host_cell_global_ids =
            Kokkos::create_mirror_view( poly_list.cell_global_ids );
        Kokkos::deep_copy( host_cell_global_ids, poly_list.cell_global_ids );
        auto host_adjacent_cells =
            Kokkos::create_mirror_view( poly_list.adjacent_cells );
        Kokkos::deep_copy( host_adjacent_cells, poly_list.adjacent_cells );
        auto host_adjacencies_per_cell =
            Kokkos::create_mirror_view( poly_list.adjacencies_per_cell );
        Kokkos::deep_copy( host_adjacencies_per_cell,
                           poly_list.adjacencies_per_cell );
        for ( unsigned i = 0; i < SIZE_1; ++i )
        {
            TEST_EQUALITY( host_cell_global_ids( i ), i + OFFSET );
            TEST_EQUALITY( host_adjacent_cells( i ), i );
            TEST_EQUALITY( host_adjacencies_per_cell( i ), 1 );
        }
    }
}

template <class UserApplication>
void test_single_topology_dof( UserApplication &user_app,
                               Teuchos::FancyOStream &out, bool &success )
{
    // Create a map.
    std::string discretization_type;
    auto dof_map = user_app.getDOFMap( discretization_type );

    // Check the map.
    TEST_EQUALITY( dof_map.object_dof_ids.rank(), 2 );
    auto host_global_dof_ids =
        Kokkos::create_mirror_view( dof_map.global_dof_ids );
    Kokkos::deep_copy( host_global_dof_ids, dof_map.global_dof_ids );
    auto host_object_dof_ids =
        Kokkos::create_mirror_view( dof_map.object_dof_ids );
    Kokkos::deep_copy( host_object_dof_ids, dof_map.object_dof_ids );
    for ( unsigned i = 0; i < SIZE_1; ++i )
    {
        TEST_EQUALITY( host_global_dof_ids( i ), i + OFFSET );
        for ( unsigned d = 0; d < SIZE_2; ++d )
            TEST_EQUALITY( host_object_dof_ids( i, d ), i + d + OFFSET );
    }
    TEST_EQUALITY( discretization_type, "unit_test_discretization" );
}

template <class UserApplication>
void test_multiple_topology_dof( UserApplication &user_app,
                                 Teuchos::FancyOStream &out, bool &success )
{
    // Create a map.
    std::string discretization_type;
    auto dof_map = user_app.getDOFMap( discretization_type );

    // Check the map.
    TEST_EQUALITY( dof_map.object_dof_ids.rank(), 1 );
    auto host_global_dof_ids =
        Kokkos::create_mirror_view( dof_map.global_dof_ids );
    Kokkos::deep_copy( host_global_dof_ids, dof_map.global_dof_ids );
    auto host_object_dof_ids =
        Kokkos::create_mirror_view( dof_map.object_dof_ids );
    Kokkos::deep_copy( host_object_dof_ids, dof_map.object_dof_ids );
    auto host_dofs_per_object =
        Kokkos::create_mirror_view( dof_map.dofs_per_object );
    Kokkos::deep_copy( host_dofs_per_object, dof_map.dofs_per_object );
    for ( unsigned i = 0; i < SIZE_1; ++i )
    {
        TEST_EQUALITY( host_global_dof_ids( i ), i + OFFSET );
        TEST_EQUALITY( host_object_dof_ids( i ), i + OFFSET );
        TEST_EQUALITY( host_dofs_per_object( i ), SIZE_2 );
    }
    TEST_EQUALITY( discretization_type, "unit_test_discretization" );
}

template <class UserApplication>
void test_field_push_pull( UserApplication &user_app,
                           Teuchos::FancyOStream &out, bool &success )
{
    using ExecutionSpace = typename UserApplication::ExecutionSpace;

    // Create a field.
    auto field_1 = user_app.getField( FIELD_NAME );

    // Put some data in the field.
    auto fill_field = KOKKOS_LAMBDA( const size_t i )
    {
        for ( unsigned d = 0; d < SPACE_DIM; ++d )
            field_1.dofs( i, d ) = i + d;
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, SIZE_1 ),
                          fill_field );
    Kokkos::fence();

    // Push the field into the app.
    user_app.pushField( FIELD_NAME, field_1 );

    // Create a second field.
    auto field_2 = user_app.getField( FIELD_NAME );

    // Pull the field out of the app.
    user_app.pullField( FIELD_NAME, field_2 );

    // Check the pulled field.
    auto host_dofs = Kokkos::create_mirror_view( field_2.dofs );
    Kokkos::deep_copy( host_dofs, field_2.dofs );
    for ( unsigned i = 0; i < SIZE_1; ++i )
        for ( unsigned d = 0; d < SPACE_DIM; ++d )
            TEST_EQUALITY( host_dofs( i, d ), i + d );
}

template <class UserApplication>
void test_field_eval( UserApplication &user_app, Teuchos::FancyOStream &out,
                      bool &success )
{
    using ExecutionSpace = typename UserApplication::ExecutionSpace;

    // Create an evaluation set.
    auto eval_set = DataTransferKit::InputAllocators<
        Kokkos::LayoutLeft, ExecutionSpace>::allocateEvaluationSet( SIZE_1,
                                                                    SPACE_DIM );
    auto fill_eval_set = KOKKOS_LAMBDA( const size_t i )
    {
        for ( unsigned d = 0; d < SPACE_DIM; ++d )
            eval_set.evaluation_points( i, d ) = i + d;
        eval_set.object_ids( i ) = i;
    };
    Kokkos::parallel_for( Kokkos::RangePolicy<ExecutionSpace>( 0, SIZE_1 ),
                          fill_eval_set );
    Kokkos::fence();

    // Create a field.
    auto field = user_app.getField( FIELD_NAME );

    // Evaluate the field.
    user_app.evaluateField( FIELD_NAME, eval_set, field );

    // Check the evaluation.
    auto host_dofs = Kokkos::create_mirror_view( field.dofs );
    Kokkos::deep_copy( host_dofs, field.dofs );
    auto host_points = Kokkos::create_mirror_view( eval_set.evaluation_points );
    Kokkos::deep_copy( host_points, eval_set.evaluation_points );
    auto host_object_ids = Kokkos::create_mirror_view( eval_set.object_ids );
    Kokkos::deep_copy( host_object_ids, eval_set.object_ids );
    for ( unsigned i = 0; i < SIZE_1; ++i )
        for ( unsigned d = 0; d < SPACE_DIM; ++d )
            TEST_EQUALITY( host_dofs( i, d ),
                           host_points( i, d ) + host_object_ids( i ) );
}

template <class UserApplication>
void test_missing_function( UserApplication &user_app,
                            Teuchos::FancyOStream &out, bool &success )
{
    // Get a node list. This should throw because the function is missing.
    bool caught_exception = false;
    try
    {
        auto node_list = user_app.getNodeList();
    }
    catch ( DataTransferKit::DataTransferKitException &e )
    {
        caught_exception = true;
    }
    TEST_ASSERT( caught_exception );
}

template <class UserApplication>
void test_too_many_functions( UserApplication &user_app,
                              Teuchos::FancyOStream &out, bool &success )
{
    // Get a dof id map. We registered both mixed and single topology
    // function so this will fail.
    bool caught_exception = false;
    try
    {
        std::string discretization_type;
        auto dof_map = user_app.getDOFMap( discretization_type );
    }
    catch ( DataTransferKit::DataTransferKitException &e )
    {
        caught_exception = true;
    }
    TEST_ASSERT( caught_exception );
}

#endif // DTK_TESTAPPLICATIONHELPERS_HPP
