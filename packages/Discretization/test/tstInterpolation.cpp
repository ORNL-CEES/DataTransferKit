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

#include "MeshGenerator.hpp"
#include <DTK_Interpolation.hpp>
#include <DTK_Mesh.hpp>
#include <DTK_Types.h>

#include <Teuchos_UnitTestHarness.hpp>

template <typename DeviceType>
Kokkos::View<DataTransferKit::Coordinate *[3], DeviceType>
getPointsCoord3D( MPI_Comm comm ) {
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    // Create the points we want to search
    unsigned int const n_points = comm_rank < 2 ? 5 : 0;
    Kokkos::View<DataTransferKit::Coordinate * [3], DeviceType> points_coord(
        "points_coord", n_points );
    auto points_coord_host = Kokkos::create_mirror_view( points_coord );
    if ( comm_rank == 0 )
    {
        // Point in the center of the cell
        points_coord_host( 0, 0 ) = 0.5;
        points_coord_host( 0, 1 ) = 0.5;
        points_coord_host( 0, 2 ) = 0.5;
        // Point off the center of the cell
        points_coord_host( 1, 0 ) = 1.25;
        points_coord_host( 1, 1 ) = 2.75;
        points_coord_host( 1, 2 ) = 3.25;
        // Point on a face
        points_coord_host( 2, 0 ) = 2.75;
        points_coord_host( 2, 1 ) = 2;
        points_coord_host( 2, 2 ) = 3.25;
        // Point on an edge
        points_coord_host( 3, 0 ) = 2.5;
        points_coord_host( 3, 1 ) = 2;
        points_coord_host( 3, 2 ) = 3;
        // Point on a vertex
        points_coord_host( 4, 0 ) = 2;
        points_coord_host( 4, 1 ) = 2;
        points_coord_host( 4, 2 ) = 2;
    }
    if ( comm_rank == 1 )
    {
        // Point in the center of the cell
        points_coord_host( 0, 0 ) = 1.5;
        points_coord_host( 0, 1 ) = 1.5;
        points_coord_host( 0, 2 ) = 1.5;
        // Point off the center of the cell
        points_coord_host( 1, 0 ) = 2.25;
        points_coord_host( 1, 1 ) = 3.75;
        points_coord_host( 1, 2 ) = 4.25;
        // Point on a face
        points_coord_host( 2, 0 ) = 3.75;
        points_coord_host( 2, 1 ) = 3;
        points_coord_host( 2, 2 ) = 4.25;
        // Point on an edge
        points_coord_host( 3, 0 ) = 3.5;
        points_coord_host( 3, 1 ) = 3;
        points_coord_host( 3, 2 ) = 4;
        // Point on a vertex
        points_coord_host( 4, 0 ) = 3;
        points_coord_host( 4, 1 ) = 3;
        points_coord_host( 4, 2 ) = 3;
    }
    Kokkos::deep_copy( points_coord, points_coord_host );

    return points_coord;
}

template <typename DeviceType>
Kokkos::View<DataTransferKit::Coordinate *[3], DeviceType> getPointsCoord3DHdiv(
    MPI_Comm comm ) {
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    // Create the points we want to search
    unsigned int const n_points = comm_rank < 2 ? 5 : 0;
    Kokkos::View<DataTransferKit::Coordinate * [3], DeviceType> points_coord(
        "points_coord", n_points );
    auto points_coord_host = Kokkos::create_mirror_view( points_coord );
    if ( comm_rank == 0 )
    {
        // Point in the center of the cell
        points_coord_host( 0, 0 ) = 0.5;
        points_coord_host( 0, 1 ) = 0.5;
        points_coord_host( 0, 2 ) = 0.5;
        // Point off the center of the cell
        points_coord_host( 1, 0 ) = 1.25;
        points_coord_host( 1, 1 ) = 2.75;
        points_coord_host( 1, 2 ) = 3.25;
        // Point off the center of the cell
        points_coord_host( 2, 0 ) = 2.75;
        points_coord_host( 2, 1 ) = 2.1;
        points_coord_host( 2, 2 ) = 3.25;
        // Point off the center of the cell
        points_coord_host( 3, 0 ) = 2.5;
        points_coord_host( 3, 1 ) = 2.2;
        points_coord_host( 3, 2 ) = 3.3;
        // Point off the center of the cell
        points_coord_host( 4, 0 ) = 2.2;
        points_coord_host( 4, 1 ) = 2.1;
        points_coord_host( 4, 2 ) = 2.3;
    }
    if ( comm_rank == 1 )
    {
        // Point in the center of the cell
        points_coord_host( 0, 0 ) = 1.5;
        points_coord_host( 0, 1 ) = 1.5;
        points_coord_host( 0, 2 ) = 1.5;
        // Point off the center of the cell
        points_coord_host( 1, 0 ) = 2.25;
        points_coord_host( 1, 1 ) = 3.75;
        points_coord_host( 1, 2 ) = 4.25;
        // Point off the center of the cell
        points_coord_host( 2, 0 ) = 3.75;
        points_coord_host( 2, 1 ) = 3.1;
        points_coord_host( 2, 2 ) = 4.25;
        // Point off the center of the cell
        points_coord_host( 3, 0 ) = 3.5;
        points_coord_host( 3, 1 ) = 3.2;
        points_coord_host( 3, 2 ) = 4.3;
        // Point off the center of the cell
        points_coord_host( 4, 0 ) = 3.3;
        points_coord_host( 4, 1 ) = 3.4;
        points_coord_host( 4, 2 ) = 3.5;
    }
    Kokkos::deep_copy( points_coord, points_coord_host );

    return points_coord;
}

template <typename DeviceType>
Kokkos::View<DataTransferKit::Coordinate *[2], DeviceType> getPointsCoord2D(
    MPI_Comm comm ) {
    int comm_size;
    MPI_Comm_size( comm, &comm_size );
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );

    // Create the points, we are looking for
    unsigned int const query_offset = 3 * ( ( comm_rank + 1 ) % comm_size );
    unsigned int constexpr n_points = 4;
    Kokkos::View<DataTransferKit::Coordinate * [2], DeviceType> points_coord(
        "points_coord", n_points );
    auto points_coord_host = Kokkos::create_mirror_view( points_coord );
    // First point
    points_coord_host( 0, 0 ) = query_offset + 1.5;
    points_coord_host( 0, 1 ) = 0.;
    // Second point
    points_coord_host( 1, 0 ) = query_offset + 1.;
    points_coord_host( 1, 1 ) = 1.5;
    // Third point
    points_coord_host( 2, 0 ) = query_offset + 1.5;
    points_coord_host( 2, 1 ) = 1.5;
    // Fourth point
    points_coord_host( 3, 0 ) = query_offset + 2.;
    points_coord_host( 3, 1 ) = 1.;
    Kokkos::deep_copy( points_coord, points_coord_host );

    return points_coord;
}

// The `out` and `success` parameters come from the Teuchos unit testing macros
// expansion.
template <int dim, typename DeviceType>
void checkReferencePoints(
    Kokkos::View<ArborX::Point *, DeviceType> phys_points,
    Kokkos::View<ArborX::Point *, DeviceType> reference_points,
    Kokkos::View<bool *, DeviceType> point_in_cell,
    std::map<std::array<double, dim>, std::vector<std::array<double, dim>>>
        &ref_sol,
    unsigned int const ref_n_points_in_cell, bool &success,
    Teuchos::FancyOStream &out )
{
    auto phys_points_host = Kokkos::create_mirror_view( phys_points );
    Kokkos::deep_copy( phys_points_host, phys_points );
    auto reference_points_host = Kokkos::create_mirror_view( reference_points );
    Kokkos::deep_copy( reference_points_host, reference_points );
    auto point_in_cell_host = Kokkos::create_mirror_view( point_in_cell );
    Kokkos::deep_copy( point_in_cell_host, point_in_cell );
    unsigned int n_points_in_cell = 0;
    for ( unsigned int i = 0; i < point_in_cell_host.extent( 0 ); ++i )
        if ( point_in_cell_host( i ) )
            ++n_points_in_cell;
    TEST_EQUALITY( n_points_in_cell, ref_n_points_in_cell );

    typedef std::array<DataTransferKit::Coordinate, dim> PtCoord;
    for ( unsigned int i = 0; i < phys_points.extent( 0 ); ++i )
    {
        if ( point_in_cell_host( i ) )
        {
            PtCoord pt_coord;
            PtCoord ref_pt;
            for ( unsigned int d = 0; d < dim; ++d )
            {
                pt_coord[d] = phys_points( i )[d];
                ref_pt[d] = reference_points_host( i )[d];
            }

            unsigned int const n_ref_points = ref_sol[pt_coord].size();
            bool pt_found;
            for ( unsigned int j = 0; j < n_ref_points; ++j )
            {
                pt_found = true;
                for ( unsigned int d = 0; d < dim; ++d )
                    if ( std::abs( ref_pt[d] - ref_sol[pt_coord][j][d] ) >
                         1e-14 )
                        pt_found = false;

                if ( pt_found == true )
                    break;
            }

            TEST_EQUALITY( pt_found, true );
        }
    }
}

// The `out` and `success` parameters come from the Teuchos unit testing macros
// expansion.
template <int dim, int ref_size, typename DeviceType>
void checkFieldValue( std::array<double, ref_size> const &ref_sol,
                      Kokkos::View<double **, DeviceType> Y, bool &success,
                      Teuchos::FancyOStream &out, double tol = 1e-14 )
{
    auto Y_host = Kokkos::create_mirror_view( Y );
    Kokkos::deep_copy( Y_host, Y );
    TEST_EQUALITY( Y.extent( 0 ), ref_size );
    unsigned int const n_fields = Y.extent( 1 );
    for ( unsigned int i = 0; i < ref_size; ++i )
        for ( unsigned int j = 0; j < n_fields; ++j )
            TEST_FLOATING_EQUALITY( ref_sol[i] + dim * j, Y_host( i, j ), tol );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Interpolation, one_topo_one_fe_three_dim,
                                   DeviceType )
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    unsigned int constexpr dim = 3;
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> coordinates;
    Kokkos::View<DataTransferKit::Coordinate * [3], DeviceType> points_coord;
    std::vector<unsigned int> n_subdivisions = {{5, 5, 3}};
    std::tie( cell_topologies, cells, coordinates ) =
        buildStructuredMesh<DeviceType>( comm, n_subdivisions );
    points_coord = getPointsCoord3D<DeviceType>( comm );
    unsigned int const n_points = points_coord.extent( 0 );

    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int const n_dofs = coordinates.extent( 0 );
    unsigned int const n_fields = 1;
    Kokkos::View<DataTransferKit::LocalOrdinal *, DeviceType> cell_dofs_ids(
        "cell_dofs_ids", cells.extent( 0 ) );

    Kokkos::parallel_for(
        "initialize_cell_dofs_ids",
        Kokkos::RangePolicy<ExecutionSpace>( 0, cells.extent( 0 ) ),
        KOKKOS_LAMBDA( int const i ) { cell_dofs_ids( i ) = cells( i ); } );
    Kokkos::fence();

    // We are now done with building the mesh and we can do the
    // interpolation
    DataTransferKit::Mesh<DeviceType> mesh( cell_topologies, cells,
                                            coordinates );
    DataTransferKit::Interpolation<DeviceType> interpolation(
        comm, mesh, points_coord, cell_dofs_ids, DTK_HGRAD );

    // We set X = x + y +z
    Kokkos::View<double **, DeviceType> X( "X", n_dofs, n_fields );
    Kokkos::parallel_for( "initialize_X",
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_dofs ),
                          KOKKOS_LAMBDA( int const i ) {
                              for ( unsigned int d = 0; d < dim; ++d )
                                  X( i, 0 ) += coordinates( i, d );
                          } );
    Kokkos::fence();

    Kokkos::View<double **, DeviceType> Y( "Y", n_points, n_fields );
    interpolation.apply( X, Y );
    if ( comm_rank == 0 )
    {
        std::array<double, 5> ref_sol = {{1.5, 7.25, 8.0, 7.5, 6.}};
        checkFieldValue<dim, 5>( ref_sol, Y, success, out );
    }
    else if ( comm_rank == 1 )
    {
        std::array<double, 5> ref_sol = {{4.5, 10.25, 11.0, 10.5, 9}};
        checkFieldValue<dim, 5>( ref_sol, Y, success, out );
    }
    else
    {
        TEST_EQUALITY( Y.extent( 0 ), 0 );
    }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Interpolation, two_topo_two_dim, DeviceType )
{
    // Test a mesh of made of Quadrilateral<4> and Triangle<3>
    // 7-----8-----9
    // |    / \    |
    // 3---4---5---6
    // |    \ /    |
    // 0-----1-----2

    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    int comm_size;
    MPI_Comm_size( comm, &comm_size );

    unsigned int constexpr dim = 2;
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> coordinates;
    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> points_coord;
    std::tie( cell_topologies, cells, coordinates ) =
        buildMixedMesh<DeviceType>( comm, 2 );
    points_coord = getPointsCoord2D<DeviceType>( comm );
    unsigned int const n_points = points_coord.extent( 0 );

    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int const n_dofs = coordinates.extent( 0 );
    unsigned int const n_fields = 2;
    Kokkos::View<DataTransferKit::LocalOrdinal *, DeviceType> cell_dofs_ids(
        "cell_dofs_ids", cells.extent( 0 ) );

    Kokkos::parallel_for(
        "initialize_cell_dofs_ids",
        Kokkos::RangePolicy<ExecutionSpace>( 0, cells.extent( 0 ) ),
        KOKKOS_LAMBDA( int const i ) { cell_dofs_ids( i ) = cells( i ); } );
    Kokkos::fence();

    // We are now done with building the mesh and we can do the interpolation
    DataTransferKit::Mesh<DeviceType> mesh( cell_topologies, cells,
                                            coordinates );
    DataTransferKit::Interpolation<DeviceType> interpolation(
        comm, mesh, points_coord, cell_dofs_ids, DTK_HGRAD );

    // We set X = x + y + 2*field_id with field_id = 0 or 1
    Kokkos::View<double **, DeviceType> X( "X", n_dofs, n_fields );
    Kokkos::parallel_for( "initialize_X",
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_dofs ),
                          KOKKOS_LAMBDA( int const i ) {
                              for ( unsigned int d = 0; d < dim; ++d )
                                  for ( unsigned int j = 0; j < n_fields; ++j )
                                  {
                                      X( i, j ) += j + coordinates( i, d );
                                  }
                          } );
    Kokkos::fence();

    Kokkos::View<double **, DeviceType> Y( "Y", n_points, n_fields );
    interpolation.apply( X, Y );

    unsigned int const query_offset = 3 * ( ( comm_rank + 1 ) % comm_size );
    std::array<double, 4> ref_sol = {{query_offset + 1.5, query_offset + 2.5,
                                      query_offset + 3., query_offset + 3.}};
    checkFieldValue<dim, 4>( ref_sol, Y, success, out );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Interpolation,
                                   one_topo_one_fe_three_dim_hdiv, DeviceType )
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> coordinates;
    Kokkos::View<DataTransferKit::Coordinate * [3], DeviceType> points_coord;
    std::vector<unsigned int> n_subdivisions = {{5, 5, 3}};
    std::tie( cell_topologies, cells, coordinates ) =
        buildStructuredMesh<DeviceType>( comm, n_subdivisions );
    points_coord = getPointsCoord3DHdiv<DeviceType>( comm );
    unsigned int const n_points = points_coord.extent( 0 );

    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int const n_dofs = 1.5 * coordinates.extent( 0 );
    unsigned int const n_fields = 1;
    Kokkos::View<DataTransferKit::LocalOrdinal *, DeviceType> cell_dofs_ids(
        "cell_dofs_ids", cells.extent( 0 ) );

    Kokkos::parallel_for(
        "initialize_cell_dofs_ids",
        Kokkos::RangePolicy<ExecutionSpace>( 0, cells.extent( 0 ) ),
        KOKKOS_LAMBDA( int const i ) { cell_dofs_ids( i ) = cells( i ); } );
    Kokkos::fence();

    // We are now done with building the mesh and we can do the
    // interpolation
    DataTransferKit::Mesh<DeviceType> mesh( cell_topologies, cells,
                                            coordinates );
    DataTransferKit::Interpolation<DeviceType> interpolation(
        comm, mesh, points_coord, cell_dofs_ids, DTK_HDIV );

    // We set X = 1. I don't know what the FE looks like so this tests just
    // check that we don't crash
    Kokkos::View<double **, DeviceType> X( "X", n_dofs, n_fields );
    Kokkos::parallel_for( "initialize_X",
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_dofs ),
                          KOKKOS_LAMBDA( int const i ) { X( i, 0 ) = 1.; } );
    Kokkos::fence();

    Kokkos::View<double **, DeviceType> Y( "Y", n_points, n_fields );
    interpolation.apply( X, Y );

    if ( comm_rank <= 1 )
    {
        TEST_EQUALITY( Y.extent( 0 ), 5 );
    }
    else
    {
        TEST_EQUALITY( Y.extent( 0 ), 0 );
    }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Interpolation,
                                   one_topo_one_fe_three_dim_point_not_found,
                                   DeviceType )
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    unsigned int constexpr dim = 3;
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> coordinates;
    Kokkos::View<DataTransferKit::Coordinate * [3], DeviceType> points_coord;
    std::vector<unsigned int> n_subdivisions = {{5, 5, 3}};
    std::tie( cell_topologies, cells, coordinates ) =
        buildStructuredMesh<DeviceType>( comm, n_subdivisions );
    points_coord = getPointsCoord3D<DeviceType>( comm );
    unsigned int const n_points = points_coord.extent( 0 );

    auto points_coord_host = Kokkos::create_mirror_view( points_coord );
    Kokkos::deep_copy( points_coord_host, points_coord );
    if ( comm_rank == 0 )
    {
        points_coord_host( 1, 0 ) = -100.;
        points_coord_host( 1, 1 ) = -100.;
        points_coord_host( 1, 2 ) = -100.;
    }
    Kokkos::deep_copy( points_coord, points_coord_host );

    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int const n_dofs = coordinates.extent( 0 );
    unsigned int const n_fields = 1;
    Kokkos::View<DataTransferKit::LocalOrdinal *, DeviceType> cell_dofs_ids(
        "cell_dofs_ids", cells.extent( 0 ) );

    Kokkos::parallel_for(
        "initialize_cell_dofs_ids",
        Kokkos::RangePolicy<ExecutionSpace>( 0, cells.extent( 0 ) ),
        KOKKOS_LAMBDA( int const i ) { cell_dofs_ids( i ) = cells( i ); } );
    Kokkos::fence();

    // We are now done with building the mesh and we can do the
    // interpolation
    DataTransferKit::Mesh<DeviceType> mesh( cell_topologies, cells,
                                            coordinates );
    DataTransferKit::Interpolation<DeviceType> interpolation(
        comm, mesh, points_coord, cell_dofs_ids, DTK_HGRAD );

    // We set X = x + y +z
    Kokkos::View<double **, DeviceType> X( "X", n_dofs, n_fields );
    Kokkos::parallel_for( "initialize_X",
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_dofs ),
                          KOKKOS_LAMBDA( int const i ) {
                              for ( unsigned int d = 0; d < dim; ++d )
                                  X( i, 0 ) += coordinates( i, d );
                          } );
    Kokkos::fence();

    Kokkos::View<double **, DeviceType> Y( "Y", n_points, n_fields );
    auto query_ids = interpolation.apply( X, Y );
    if ( comm_rank == 0 )
    {
        unsigned int constexpr array_size = 5;
        std::array<double, array_size> ref_sol = {{1.5, 8.0, 7.5, 6., 0.}};
        checkFieldValue<dim, array_size>( ref_sol, Y, success, out );

        auto query_ids_host = Kokkos::create_mirror_view( query_ids );
        Kokkos::deep_copy( query_ids_host, query_ids );
        std::array<int, array_size> ref_query_id = {{0, 2, 3, 4, -1}};
        for ( unsigned int i = 0; i < array_size; ++i )
            TEST_EQUALITY( query_ids_host( i ), ref_query_id[i] );
    }
    else if ( comm_rank == 1 )
    {
        unsigned int constexpr array_size = 5;
        std::array<double, 5> ref_sol = {{4.5, 10.25, 11.0, 10.5, 9}};
        checkFieldValue<dim, 5>( ref_sol, Y, success, out );

        auto query_ids_host = Kokkos::create_mirror_view( query_ids );
        Kokkos::deep_copy( query_ids_host, query_ids );
        std::array<int, array_size> ref_query_id = {{0, 1, 2, 3, 4}};
        for ( unsigned int i = 0; i < array_size; ++i )
            TEST_EQUALITY( query_ids_host( i ), ref_query_id[i] );
    }
    else
    {
        TEST_EQUALITY( Y.extent( 0 ), 0 );
    }
}

// Include the test macros.
#include "DataTransferKit_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    using DeviceType##NODE = typename NODE::device_type;                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                      \
        Interpolation, one_topo_one_fe_three_dim, DeviceType##NODE )           \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Interpolation, two_topo_two_dim,     \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                      \
        Interpolation, one_topo_one_fe_three_dim_hdiv, DeviceType##NODE )      \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                      \
        Interpolation, one_topo_one_fe_three_dim_point_not_found,              \
        DeviceType##NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
