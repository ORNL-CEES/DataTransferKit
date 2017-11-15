/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#include <DTK_Interpolation.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_UnitTestHarness.hpp>

template <typename DeviceType>
std::tuple<Kokkos::View<DTK_CellTopology *, DeviceType>,
           Kokkos::View<unsigned int *, DeviceType>,
           Kokkos::View<double **, DeviceType>,
           Kokkos::View<double * [3], DeviceType>>
buildHexMesh( Teuchos::RCP<const Teuchos::Comm<int>> comm )
{
    int const comm_rank = comm->getRank();
    unsigned int constexpr dim = 3;
    unsigned int constexpr i_max = 3;
    unsigned int constexpr j_max = 5;
    unsigned int constexpr k_max = 5;
    unsigned int constexpr n_local_cells =
        ( i_max - 1 ) * ( j_max - 1 ) * ( k_max - 1 );
    unsigned int constexpr n_vertices = i_max * j_max * k_max;

    // Create the Kokkos::View of the coordinates
    Kokkos::View<double **, DeviceType> coordinates( "coordinates", n_vertices,
                                                     dim );
    auto coordinates_host = Kokkos::create_mirror_view( coordinates );
    unsigned int n = 0;
    for ( unsigned int i = 0; i < i_max; ++i )
    {
        unsigned int const m = 2 * comm_rank + i;
        for ( unsigned int j = 0; j < j_max; ++j )
            for ( unsigned int k = 0; k < k_max; ++k )
            {
                coordinates_host( n, 0 ) = static_cast<double>( k );
                coordinates_host( n, 1 ) = static_cast<double>( j );
                coordinates_host( n, 2 ) = static_cast<double>( m );
                ++n;
            }
    }
    Kokkos::deep_copy( coordinates, coordinates_host );

    // Create the Kokkos::View of the coordinates
    unsigned int constexpr n_vertices_per_cell = 8;
    Kokkos::View<unsigned int *, DeviceType> cells(
        "cells", n_local_cells * n_vertices_per_cell );
    auto cells_host = Kokkos::create_mirror_view( cells );
    n = 0;
    for ( unsigned int i = 0; i < i_max - 1; ++i )
        for ( unsigned int j = 0; j < j_max - 1; ++j )
            for ( unsigned int k = 0; k < k_max - 1; ++k )
            {
                cells_host( n++ ) = k + j * k_max + i * j_max * k_max;
                cells_host( n++ ) = ( k + 1 ) + j * k_max + i * j_max * k_max;
                cells_host( n++ ) =
                    ( k + 1 ) + ( j + 1 ) * k_max + i * j_max * k_max;
                cells_host( n++ ) = k + ( j + 1 ) * k_max + i * j_max * k_max;
                cells_host( n++ ) = k + j * k_max + ( i + 1 ) * j_max * k_max;
                cells_host( n++ ) =
                    ( k + 1 ) + j * k_max + ( i + 1 ) * j_max * k_max;
                cells_host( n++ ) =
                    ( k + 1 ) + ( j + 1 ) * k_max + ( i + 1 ) * j_max * k_max;
                cells_host( n++ ) =
                    k + ( j + 1 ) * k_max + ( i + 1 ) * j_max * k_max;
            }
    Kokkos::deep_copy( cells, cells_host );

    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view(
        "cell_topologies", n_local_cells );
    Kokkos::deep_copy( cell_topologies_view, DTK_HEX_8 );

    // Create the points we want to search
    unsigned int const n_points = comm_rank < 2 ? 5 : 0;
    Kokkos::View<double * [3], DeviceType> points_coord( "points_coord",
                                                         n_points );
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

    return std::make_tuple( cell_topologies_view, cells, coordinates,
                            points_coord );
}

template <typename DeviceType>
std::tuple<Kokkos::View<DTK_CellTopology *, DeviceType>,
           Kokkos::View<unsigned int *, DeviceType>,
           Kokkos::View<double **, DeviceType>,
           Kokkos::View<double **, DeviceType>>
buildMixedMesh( Teuchos::RCP<const Teuchos::Comm<int>> comm )
{
    int const comm_rank = comm->getRank();
    int const comm_size = comm->getSize();
    // Build the mesh
    unsigned int constexpr dim = 2;
    unsigned int constexpr n_local_quad_cells = 4;
    unsigned int constexpr n_local_tri_cells = 2;
    unsigned int constexpr n_local_cells =
        n_local_quad_cells + n_local_tri_cells;
    unsigned int const offset_mesh = 3 * comm_rank;
    unsigned int constexpr n_vertices = 10;

    // Create the Kokkos::View of the coordinates
    Kokkos::View<double **, DeviceType> coordinates( "coordinates", n_vertices,
                                                     dim );
    auto coordinates_host = Kokkos::create_mirror_view( coordinates );
    // Y=0 points
    coordinates_host( 0, 0 ) = offset_mesh;
    coordinates_host( 0, 1 ) = 0.;
    coordinates_host( 1, 0 ) = offset_mesh + 1.5;
    coordinates_host( 1, 1 ) = 0.;
    coordinates_host( 2, 0 ) = offset_mesh + 3.;
    coordinates_host( 2, 1 ) = 0.;
    // Y=1 points
    coordinates_host( 3, 0 ) = offset_mesh;
    coordinates_host( 3, 1 ) = 1.;
    coordinates_host( 4, 0 ) = offset_mesh + 1;
    coordinates_host( 4, 1 ) = 1.;
    coordinates_host( 5, 0 ) = offset_mesh + 2;
    coordinates_host( 5, 1 ) = 1.;
    coordinates_host( 6, 0 ) = offset_mesh + 3;
    coordinates_host( 6, 1 ) = 1.;
    // Y=2 points
    coordinates_host( 7, 0 ) = offset_mesh;
    coordinates_host( 7, 1 ) = 2.;
    coordinates_host( 8, 0 ) = offset_mesh + 1.5;
    coordinates_host( 8, 1 ) = 2.;
    coordinates_host( 9, 0 ) = offset_mesh + 3.;
    coordinates_host( 9, 1 ) = 2.;
    Kokkos::deep_copy( coordinates, coordinates_host );

    // Create the Kokkos::View of the coordinates
    Kokkos::View<unsigned int *, DeviceType> cells(
        "cells", 4 * n_local_quad_cells + 3 * n_local_tri_cells );
    auto cells_host = Kokkos::create_mirror_view( cells );
    unsigned int n = 0;
    // First cell
    cells_host( n++ ) = 0;
    cells_host( n++ ) = 1;
    cells_host( n++ ) = 4;
    cells_host( n++ ) = 3;
    // Second cell
    cells_host( n++ ) = 1;
    cells_host( n++ ) = 5;
    cells_host( n++ ) = 4;
    // Third cell
    cells_host( n++ ) = 1;
    cells_host( n++ ) = 2;
    cells_host( n++ ) = 6;
    cells_host( n++ ) = 5;
    // Fourth cell
    cells_host( n++ ) = 3;
    cells_host( n++ ) = 4;
    cells_host( n++ ) = 8;
    cells_host( n++ ) = 7;
    // Fifth cell
    cells_host( n++ ) = 5;
    cells_host( n++ ) = 8;
    cells_host( n++ ) = 4;
    // Sixth cell
    cells_host( n++ ) = 5;
    cells_host( n++ ) = 6;
    cells_host( n++ ) = 9;
    cells_host( n++ ) = 8;
    Kokkos::deep_copy( cells, cells_host );

    // Create the Kokkos::View of the topologies
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view(
        "cell_topologies", n_local_cells );
    auto cell_topologies_view_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    cell_topologies_view_host( 0 ) = DTK_QUAD_4;
    cell_topologies_view_host( 1 ) = DTK_TRI_3;
    cell_topologies_view_host( 2 ) = DTK_QUAD_4;
    cell_topologies_view_host( 3 ) = DTK_QUAD_4;
    cell_topologies_view_host( 4 ) = DTK_TRI_3;
    cell_topologies_view_host( 5 ) = DTK_QUAD_4;
    Kokkos::deep_copy( cell_topologies_view, cell_topologies_view_host );

    // Create the points, we are looking for
    unsigned int const query_offset = 3 * ( ( comm_rank + 1 ) % comm_size );
    unsigned int constexpr n_points = 4;
    Kokkos::View<double **, DeviceType> points_coord( "points_coord", n_points,
                                                      2 );
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
    Kokkos::deep_copy( points_coord_host, points_coord );

    return std::make_tuple( cell_topologies_view, cells, coordinates,
                            points_coord );
}

// The `out` and `success` parameters come from the Teuchos unit testing macros
// expansion.
template <int dim, typename DeviceType>
void checkReferencePoints(
    Kokkos::View<DataTransferKit::Point *, DeviceType> phys_points,
    Kokkos::View<DataTransferKit::Point *, DeviceType> reference_points,
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

    typedef std::array<double, dim> PtCoord;
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
                      Teuchos::FancyOStream &out )
{
    auto Y_host = Kokkos::create_mirror_view( Y );
    Kokkos::deep_copy( Y_host, Y );
    TEST_EQUALITY( Y.extent( 0 ), ref_size );
    unsigned int const n_fields = Y.extent( 1 );
    for ( unsigned int i = 0; i < ref_size; ++i )
        for ( unsigned int j = 0; j < n_fields; ++j )
            TEST_FLOATING_EQUALITY( ref_sol[i] + dim * j, Y_host( i, j ),
                                    1e-14 );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Interpolation, one_topo_three_dim,
                                   DeviceType )
{
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int const comm_rank = comm->getRank();
    unsigned int constexpr dim = 3;
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<double **, DeviceType> coordinates;
    Kokkos::View<double * [3], DeviceType> points_coord;
    std::tie( cell_topologies_view, cells, coordinates, points_coord ) =
        buildHexMesh<DeviceType>( comm );

    // We are now done with building the mesh and we can do the interpolation
    DataTransferKit::Interpolation<DeviceType> interpolation(
        comm, cell_topologies_view, cells, coordinates, points_coord,
        DTK_HGRAD );

    Kokkos::View<DataTransferKit::Point *, DeviceType> phys_points(
        "phys_points" );
    Kokkos::View<DataTransferKit::Point *, DeviceType> reference_points(
        "reference_points" );
    Kokkos::View<bool *, DeviceType> point_in_cell( "pt_in_cell" );
    interpolation.getReferencePoints( phys_points, reference_points,
                                      point_in_cell );

    // Check the number of points found on each processor
    if ( comm_rank == 0 )
    {
        TEST_EQUALITY( reference_points.extent( 0 ), 16 );
    }
    else if ( comm_rank == 1 )
    {
        TEST_EQUALITY( reference_points.extent( 0 ), 16 );
    }
    else
    {
        TEST_EQUALITY( reference_points.extent( 0 ), 0 );
    }

    // Reference solution
    typedef std::array<double, dim> PtCoord;
    std::map<PtCoord, std::vector<PtCoord>> ref_sol;
    // First point
    PtCoord pt_1 = {{0.5, 0.5, 0.5}};
    std::vector<PtCoord> ref_frame_1 = {{{0., 0., 0.}}};
    ref_sol[pt_1] = ref_frame_1;
    // Second point
    PtCoord pt_2 = {{1.25, 2.75, 3.25}};
    std::vector<PtCoord> ref_frame_2 = {{{-0.5, 0.5, -0.5}}};
    ref_sol[pt_2] = ref_frame_2;
    // Third point
    PtCoord pt_3 = {{2.75, 2., 3.25}};
    std::vector<PtCoord> ref_frame_3 = {{{0.5, -1, -0.5}}, {{0.5, 1, -0.5}}};
    ref_sol[pt_3] = ref_frame_3;
    // Fourth point
    PtCoord pt_4 = {{2.5, 2., 3.}};
    std::vector<PtCoord> ref_frame_4 = {
        {{0., -1., -1.}}, {{0., 1., -1.}}, {{0., 1., 1.}}, {{0., -1., 1}}};
    ref_sol[pt_4] = ref_frame_4;
    // Fifth point
    PtCoord pt_5 = {{2., 2., 2.}};
    std::vector<PtCoord> ref_frame_5 = {
        {{-1., -1., 1.}},  {{-1., 1., 1}},  {{1., -1., 1.}}, {{1., 1., 1.}},
        {{-1., -1., -1.}}, {{-1., 1, -1.}}, {{1, -1., -1.}}, {{1., 1., -1.}}};
    ref_sol[pt_5] = ref_frame_5;
    // Sixth point
    PtCoord pt_6 = {{1.5, 1.5, 1.5}};
    std::vector<PtCoord> ref_frame_6 = {{{0., 0., 0.}}};
    ref_sol[pt_6] = ref_frame_6;
    // Seventh point
    PtCoord pt_7 = {{2.25, 3.75, 4.25}};
    std::vector<PtCoord> ref_frame_7 = {{{-0.5, 0.5, -0.5}}};
    ref_sol[pt_7] = ref_frame_7;
    // Eigth point
    PtCoord pt_8 = {{3.75, 3., 4.25}};
    std::vector<PtCoord> ref_frame_8 = {{{0.5, -1., -0.5}}, {{0.5, 1., -0.5}}};
    ref_sol[pt_8] = ref_frame_8;
    // Nineth point
    PtCoord pt_9 = {{3.5, 3., 4.}};
    std::vector<PtCoord> ref_frame_9 = {
        {{0., -1., 1.}}, {{0., 1., 1.}}, {{0., -1., -1.}}, {{0., 1., -1.}}};
    ref_sol[pt_9] = ref_frame_9;
    // Tenth point
    PtCoord pt_10 = {{3., 3., 3.}};
    std::vector<PtCoord> ref_frame_10 = {
        {{-1., -1., -1.}}, {{-1., -1., 1.}}, {{-1., 1., 1.}},  {{1., 1., 1.}},
        {{1., -1., -1.}},  {{1., 1., -1.}},  {{-1., 1., -1.}}, {{1., -1., 1.}}};
    ref_sol[pt_10] = ref_frame_10;

    auto phys_points_host = Kokkos::create_mirror_view( phys_points );
    Kokkos::deep_copy( phys_points_host, phys_points );
    auto reference_points_host = Kokkos::create_mirror_view( reference_points );
    Kokkos::deep_copy( reference_points_host, reference_points );
    auto point_in_cell_host = Kokkos::create_mirror_view( point_in_cell );
    Kokkos::deep_copy( point_in_cell_host, point_in_cell );
    unsigned int ref_n_points_in_cell = 0;
    if ( comm_rank < 2 )
        ref_n_points_in_cell = 16;
    else
        ref_n_points_in_cell = 0;

    // Check the results
    checkReferencePoints<dim>( phys_points, reference_points, point_in_cell,
                               ref_sol, ref_n_points_in_cell, success, out );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Interpolation, two_topo_two_dim, DeviceType )
{
    // Test a mesh of made of Quadrilateral<4> and Triangle<3>
    // 7-----8-----9
    // |    / \    |
    // 3---4---5---6
    // |    \ /    |
    // 0-----1-----2
    //
    // The pattern above is repeated on each processors
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int const comm_rank = comm->getRank();
    int const comm_size = comm->getSize();
    unsigned int constexpr dim = 2;
    unsigned int const query_offset = 3 * ( ( comm_rank + 1 ) % comm_size );
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<double **, DeviceType> coordinates;
    Kokkos::View<double **, DeviceType> points_coord;
    std::tie( cell_topologies_view, cells, coordinates, points_coord ) =
        buildMixedMesh<DeviceType>( comm );

    // We are now done with building the mesh and we can do the interpolation
    DataTransferKit::Interpolation<DeviceType> interpolation(
        comm, cell_topologies_view, cells, coordinates, points_coord,
        DTK_HGRAD );

    Kokkos::View<DataTransferKit::Point *, DeviceType> phys_points(
        "phys_points" );
    Kokkos::View<DataTransferKit::Point *, DeviceType> reference_points(
        "reference_points" );
    Kokkos::View<bool *, DeviceType> point_in_cell( "pt_in_cell" );
    interpolation.getReferencePoints( phys_points, reference_points,
                                      point_in_cell );

    // Build the reference solution
    unsigned int const ref_n_points_in_cell = 9;
    typedef std::array<double, dim> PtCoord;
    std::map<PtCoord, std::vector<PtCoord>> ref_sol;
    // First point
    PtCoord pt_1 = {{query_offset + 1.5, 0.}};
    std::vector<PtCoord> ref_frame_1 = {{{1., -1.}}, {{-1., -1.}}, {{0., 0.}}};
    ref_sol[pt_1] = ref_frame_1;
    // Second point
    PtCoord pt_2 = {{query_offset + 1., 1.5}};
    std::vector<PtCoord> ref_frame_2 = {{{0.6, 0.}}};
    ref_sol[pt_2] = ref_frame_2;
    // Third point
    PtCoord pt_3 = {{query_offset + 2., 1.}};
    std::vector<PtCoord> ref_frame_3 = {
        {{-1., 1.}}, {{-1., -1.}}, {{0., 0.}}, {{1., 0.}}};
    ref_sol[pt_3] = ref_frame_3;
    // Fourth point
    PtCoord pt_4 = {{query_offset + 1.5, 1.5}};
    std::vector<PtCoord> ref_frame_4 = {{{0.5, 0.25}}};
    ref_sol[pt_4] = ref_frame_4;

    // Check the results
    checkReferencePoints<dim>( phys_points, reference_points, point_in_cell,
                               ref_sol, ref_n_points_in_cell, success, out );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Interpolation,
                                   one_topo_three_dim_interpolation,
                                   DeviceType )
{
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int const comm_rank = comm->getRank();
    unsigned int constexpr dim = 3;
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<double **, DeviceType> coordinates;
    Kokkos::View<double * [dim], DeviceType> points_coord;
    std::tie( cell_topologies_view, cells, coordinates, points_coord ) =
        buildHexMesh<DeviceType>( comm );
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

    // We are now done with building the mesh and we can do the interpolation
    DataTransferKit::Interpolation<DeviceType> interpolation(
        comm, cell_topologies_view, cells, coordinates, points_coord,
        cell_dofs_ids, DTK_HGRAD );

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
        std::array<double, 0> ref_sol;
        checkFieldValue<dim, 0>( ref_sol, Y, success, out );
    }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Interpolation,
                                   two_topo_two_dim_interpolation, DeviceType )
{
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int const comm_rank = comm->getRank();
    int const comm_size = comm->getSize();
    unsigned int constexpr dim = 2;
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<double **, DeviceType> coordinates;
    Kokkos::View<double **, DeviceType> points_coord;
    std::tie( cell_topologies_view, cells, coordinates, points_coord ) =
        buildMixedMesh<DeviceType>( comm );
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
    DataTransferKit::Interpolation<DeviceType> interpolation(
        comm, cell_topologies_view, cells, coordinates, points_coord,
        cell_dofs_ids, DTK_HGRAD );

    // We set X = x + y + 2*field_id with field_id = 0 or 1
    Kokkos::View<double **, DeviceType> X( "X", n_dofs, n_fields );
    Kokkos::parallel_for( "initialize_X",
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_dofs ),
                          KOKKOS_LAMBDA( int const i ) {
                              for ( unsigned int d = 0; d < dim; ++d )
                                  for ( unsigned int j = 0; j < n_fields; ++j )
                                      X( i, j ) += j + coordinates( i, d );
                          } );
    Kokkos::fence();

    Kokkos::View<double **, DeviceType> Y( "Y", n_points, n_fields );
    interpolation.apply( X, Y );

    unsigned int const query_offset = 3 * ( ( comm_rank + 1 ) % comm_size );
    std::array<double, 4> ref_sol = {{query_offset + 1.5, query_offset + 2.5,
                                      query_offset + 3., query_offset + 3.}};
    checkFieldValue<dim, 4>( ref_sol, Y, success, out );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Interpolation, check_throw, DeviceType )
{
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<double **, DeviceType> coordinates;
    Kokkos::View<double **, DeviceType> points_coord;
    Kokkos::View<DataTransferKit::LocalOrdinal *, DeviceType> cell_dofs_ids;
    Kokkos::View<unsigned int *, DeviceType> fe_order_view;
    Kokkos::View<DTK_Quadrature *, DeviceType> quadrature_view;

    TEST_THROW( DataTransferKit::Interpolation<DeviceType> interpolation(
                    comm, cell_topologies_view, cells, coordinates,
                    points_coord, DTK_HDIV ),
                DataTransferKit::DataTransferKitException );

    TEST_THROW( DataTransferKit::Interpolation<DeviceType> interpolation(
                    comm, cell_topologies_view, cells, coordinates,
                    points_coord, cell_dofs_ids, DTK_HDIV ),
                DataTransferKit::DataTransferKitException );

    TEST_THROW( DataTransferKit::Interpolation<DeviceType> interpolation(
                    comm, cell_topologies_view, cells, coordinates,
                    points_coord, fe_order_view, quadrature_view, DTK_HGRAD ),
                DataTransferKit::DataTransferKitException );

    TEST_THROW( DataTransferKit::Interpolation<DeviceType> interpolation(
                    comm, cell_topologies_view, cells, coordinates,
                    points_coord, cell_dofs_ids, fe_order_view, quadrature_view,
                    DTK_HGRAD ),
                DataTransferKit::DataTransferKitException );
}

// Include the test macros.
#include "DataTransferKitDiscretization_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    using DeviceType##NODE = typename NODE::device_type;                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Interpolation, one_topo_three_dim,   \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Interpolation, two_topo_two_dim,     \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                      \
        Interpolation, one_topo_three_dim_interpolation, DeviceType##NODE )    \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                      \
        Interpolation, two_topo_two_dim_interpolation, DeviceType##NODE )      \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Interpolation, check_throw,          \
                                          DeviceType##NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
