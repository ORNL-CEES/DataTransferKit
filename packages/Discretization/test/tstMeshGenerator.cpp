/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#include "MeshGenerator.hpp"
#include <Teuchos_UnitTestHarness.hpp>

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MeshGenerator, structured, DeviceType )
{
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();

    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<double **, DeviceType> coordinates;

    // 2D test
    std::vector<unsigned int> n_subdivisions = {{4, 3}};
    unsigned int n_cells = 1;
    for ( auto n_sub : n_subdivisions )
        n_cells *= n_sub;
    unsigned int constexpr dim_2 = 2;
    unsigned int constexpr n_vertices_per_quad = 4;

    int const comm_rank = comm->getRank();
    double offset = n_subdivisions[1] * comm_rank;
    std::vector<std::array<std::array<double, dim_2>, n_vertices_per_quad>>
        quad_mesh_ref( n_cells );
    unsigned int n = 0;
    for ( unsigned int i = 0; i < n_subdivisions[1]; ++i )
        for ( unsigned int j = 0; j < n_subdivisions[0]; ++j )
        {
            quad_mesh_ref[n][0][0] = j;
            quad_mesh_ref[n][0][1] = i + offset;

            quad_mesh_ref[n][1][0] = j + 1;
            quad_mesh_ref[n][1][1] = i + offset;

            quad_mesh_ref[n][2][0] = j + 1;
            quad_mesh_ref[n][2][1] = i + 1 + offset;

            quad_mesh_ref[n][3][0] = j;
            quad_mesh_ref[n][3][1] = i + 1 + offset;

            ++n;
        }

    std::tie( cell_topologies_view, cells, coordinates ) =
        buildStructuredMesh<DeviceType>( comm, n_subdivisions );
    TEST_EQUALITY( cell_topologies_view.extent( 0 ), n_cells );
    TEST_EQUALITY( cells.extent( 0 ), n_cells * n_vertices_per_quad );

    auto cell_topologies_view_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    Kokkos::deep_copy( cell_topologies_view_host, cell_topologies_view );
    for ( unsigned int i = 0; i < n_cells; ++i )
        TEST_EQUALITY( cell_topologies_view_host( i ), DTK_QUAD_4 );

    auto cells_host = Kokkos::create_mirror_view( cells );
    Kokkos::deep_copy( cells_host, cells );
    auto coordinates_host = Kokkos::create_mirror_view( coordinates );
    Kokkos::deep_copy( coordinates_host, coordinates );
    n = 0;
    for ( unsigned int i = 0; i < n_cells; ++i )
    {
        for ( unsigned int j = 0; j < n_vertices_per_quad; ++j )
        {
            for ( unsigned int k = 0; k < dim_2; ++k )
            {
                unsigned int coord_pos = cells_host( n );
                TEST_EQUALITY( coordinates_host( coord_pos, k ),
                               quad_mesh_ref[i][j][k] );
            }

            ++n;
        }
    }

    // 3D test
    n_subdivisions = {{2, 3, 4}};
    n_cells = 1;
    for ( auto n_sub : n_subdivisions )
        n_cells *= n_sub;
    unsigned int constexpr dim_3 = 3;
    unsigned int constexpr n_vertices_per_hex = 8;

    offset = n_subdivisions[2] * comm_rank;
    std::vector<std::array<std::array<double, dim_3>, n_vertices_per_hex>>
        hex_mesh_ref( n_cells );
    n = 0;
    for ( unsigned int i = 0; i < n_subdivisions[2]; ++i )
        for ( unsigned int j = 0; j < n_subdivisions[1]; ++j )
            for ( unsigned int k = 0; k < n_subdivisions[0]; ++k )
            {
                hex_mesh_ref[n][0][0] = k;
                hex_mesh_ref[n][0][1] = j;
                hex_mesh_ref[n][0][2] = i + offset;

                hex_mesh_ref[n][1][0] = k + 1;
                hex_mesh_ref[n][1][1] = j;
                hex_mesh_ref[n][1][2] = i + offset;

                hex_mesh_ref[n][2][0] = k + 1;
                hex_mesh_ref[n][2][1] = j + 1;
                hex_mesh_ref[n][2][2] = i + offset;

                hex_mesh_ref[n][3][0] = k;
                hex_mesh_ref[n][3][1] = j + 1;
                hex_mesh_ref[n][3][2] = i + offset;

                hex_mesh_ref[n][4][0] = k;
                hex_mesh_ref[n][4][1] = j;
                hex_mesh_ref[n][4][2] = i + 1 + offset;

                hex_mesh_ref[n][5][0] = k + 1;
                hex_mesh_ref[n][5][1] = j;
                hex_mesh_ref[n][5][2] = i + 1 + offset;

                hex_mesh_ref[n][6][0] = k + 1;
                hex_mesh_ref[n][6][1] = j + 1;
                hex_mesh_ref[n][6][2] = i + 1 + offset;

                hex_mesh_ref[n][7][0] = k;
                hex_mesh_ref[n][7][1] = j + 1;
                hex_mesh_ref[n][7][2] = i + 1 + offset;

                ++n;
            }

    std::tie( cell_topologies_view, cells, coordinates ) =
        buildStructuredMesh<DeviceType>( comm, n_subdivisions );
    TEST_EQUALITY( cell_topologies_view.extent( 0 ), n_cells );
    TEST_EQUALITY( cells.extent( 0 ), n_cells * n_vertices_per_hex );

    cell_topologies_view_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    Kokkos::deep_copy( cell_topologies_view_host, cell_topologies_view );
    for ( unsigned int i = 0; i < n_cells; ++i )
        TEST_EQUALITY( cell_topologies_view_host( i ), DTK_HEX_8 );

    cells_host = Kokkos::create_mirror_view( cells );
    Kokkos::deep_copy( cells_host, cells );
    coordinates_host = Kokkos::create_mirror_view( coordinates );
    Kokkos::deep_copy( coordinates_host, coordinates );
    n = 0;
    for ( unsigned int i = 0; i < n_cells; ++i )
    {
        for ( unsigned int j = 0; j < n_vertices_per_hex; ++j )
        {
            for ( unsigned int k = 0; k < dim_3; ++k )
            {
                unsigned int coord_pos = cells_host( n );
                TEST_EQUALITY( coordinates_host( coord_pos, k ),
                               hex_mesh_ref[i][j][k] );
            }

            ++n;
        }
    }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MeshGenerator, mixed, DeviceType )
{
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();

    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<double **, DeviceType> coordinates;

    // 2D test
    unsigned int constexpr dim_2 = 2;
    unsigned int const n_cells = 6;
    unsigned int n_vertices = 22;
    int const comm_rank = comm->getRank();
    double offset = 3. * comm_rank;
    std::vector<std::vector<std::array<double, dim_2>>> mesh2D_ref( n_cells );
    // First cell
    mesh2D_ref[0].resize( 4 );
    mesh2D_ref[0][0][0] = offset;
    mesh2D_ref[0][0][1] = 0.;
    mesh2D_ref[0][1][0] = 1.5 + offset;
    mesh2D_ref[0][1][1] = 0.;
    mesh2D_ref[0][2][0] = 1. + offset;
    mesh2D_ref[0][2][1] = 1.;
    mesh2D_ref[0][3][0] = offset;
    mesh2D_ref[0][3][1] = 1.;
    // Second cell
    mesh2D_ref[1].resize( 3 );
    mesh2D_ref[1][0][0] = 1.5 + offset;
    mesh2D_ref[1][0][1] = 0.;
    mesh2D_ref[1][1][0] = 2. + offset;
    mesh2D_ref[1][1][1] = 1.;
    mesh2D_ref[1][2][0] = 1. + offset;
    mesh2D_ref[1][2][1] = 1.;
    // Third cell
    mesh2D_ref[2].resize( 4 );
    mesh2D_ref[2][0][0] = 1.5 + offset;
    mesh2D_ref[2][0][1] = 0.;
    mesh2D_ref[2][1][0] = 3. + offset;
    mesh2D_ref[2][1][1] = 0.;
    mesh2D_ref[2][2][0] = 3. + offset;
    mesh2D_ref[2][2][1] = 1.;
    mesh2D_ref[2][3][0] = 2. + offset;
    mesh2D_ref[2][3][1] = 1.;
    // Fourth cell
    mesh2D_ref[3].resize( 4 );
    mesh2D_ref[3][0][0] = offset;
    mesh2D_ref[3][0][1] = 1.;
    mesh2D_ref[3][1][0] = 1. + offset;
    mesh2D_ref[3][1][1] = 1.;
    mesh2D_ref[3][2][0] = 1.5 + offset;
    mesh2D_ref[3][2][1] = 2.;
    mesh2D_ref[3][3][0] = offset;
    mesh2D_ref[3][3][1] = 2.;
    // Fifth cell
    mesh2D_ref[4].resize( 3 );
    mesh2D_ref[4][0][0] = 1. + offset;
    mesh2D_ref[4][0][1] = 1.;
    mesh2D_ref[4][1][0] = 2. + offset;
    mesh2D_ref[4][1][1] = 1.;
    mesh2D_ref[4][2][0] = 1.5 + offset;
    mesh2D_ref[4][2][1] = 2.;
    // Sixth cell
    mesh2D_ref[5].resize( 4 );
    mesh2D_ref[5][0][0] = 2. + offset;
    mesh2D_ref[5][0][1] = 1.;
    mesh2D_ref[5][1][0] = 3. + offset;
    mesh2D_ref[5][1][1] = 1.;
    mesh2D_ref[5][2][0] = 3. + offset;
    mesh2D_ref[5][2][1] = 2.;
    mesh2D_ref[5][3][0] = 1.5 + offset;
    mesh2D_ref[5][3][1] = 2.;

    std::tie( cell_topologies_view, cells, coordinates ) =
        buildMixedMesh<DeviceType>( comm, dim_2 );
    TEST_EQUALITY( cell_topologies_view.extent( 0 ), n_cells );
    TEST_EQUALITY( cells.extent( 0 ), n_vertices );

    auto cell_topologies_view_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    Kokkos::deep_copy( cell_topologies_view_host, cell_topologies_view );
    std::vector<DTK_CellTopology> cell_topology_ref = {
        {DTK_QUAD_4, DTK_TRI_3, DTK_QUAD_4, DTK_QUAD_4, DTK_TRI_3, DTK_QUAD_4}};
    for ( unsigned int i = 0; i < n_cells; ++i )
        TEST_EQUALITY( cell_topologies_view_host( i ), cell_topology_ref[i] );

    auto cells_host = Kokkos::create_mirror_view( cells );
    Kokkos::deep_copy( cells_host, cells );
    auto coordinates_host = Kokkos::create_mirror_view( coordinates );
    Kokkos::deep_copy( coordinates_host, coordinates );
    unsigned int n = 0;
    for ( unsigned int i = 0; i < n_cells; ++i )
    {
        for ( unsigned int j = 0; j < mesh2D_ref[i].size(); ++j )
        {
            for ( unsigned int k = 0; k < dim_2; ++k )
            {
                unsigned int coord_pos = cells_host( n );
                TEST_EQUALITY( coordinates_host( coord_pos, k ),
                               mesh2D_ref[i][j][k] );
            }
            ++n;
        }
    }

    // 3D test
    unsigned int constexpr dim_3 = 3;
    n_vertices = 40;
    std::vector<std::vector<std::array<double, dim_3>>> mesh3D_ref( n_cells );
    // First cell
    mesh3D_ref[0].resize( 8 );
    mesh3D_ref[0][0][0] = offset;
    mesh3D_ref[0][0][1] = 0.;
    mesh3D_ref[0][0][2] = 0.;
    mesh3D_ref[0][1][0] = 1.5 + offset;
    mesh3D_ref[0][1][1] = 0.;
    mesh3D_ref[0][0][2] = 0.;
    mesh3D_ref[0][2][0] = 1. + offset;
    mesh3D_ref[0][2][1] = 1.;
    mesh3D_ref[0][0][2] = 0.;
    mesh3D_ref[0][3][0] = offset;
    mesh3D_ref[0][3][1] = 1.;
    mesh3D_ref[0][3][2] = 0.;
    mesh3D_ref[0][4][0] = offset;
    mesh3D_ref[0][4][1] = 0.;
    mesh3D_ref[0][4][2] = 1.;
    mesh3D_ref[0][5][0] = 1.5 + offset;
    mesh3D_ref[0][5][1] = 0.;
    mesh3D_ref[0][5][2] = 1.;
    mesh3D_ref[0][6][0] = 1.5 + offset;
    mesh3D_ref[0][6][1] = 1.;
    mesh3D_ref[0][6][2] = 1.;
    mesh3D_ref[0][7][0] = offset;
    mesh3D_ref[0][7][1] = 1.;
    mesh3D_ref[0][7][2] = 1.;
    // Second cell
    mesh3D_ref[1].resize( 4 );
    mesh3D_ref[1][0][0] = 1.5 + offset;
    mesh3D_ref[1][0][1] = 0.;
    mesh3D_ref[1][0][2] = 0.;
    mesh3D_ref[1][1][0] = 2. + offset;
    mesh3D_ref[1][1][1] = 1.;
    mesh3D_ref[1][1][2] = 0.;
    mesh3D_ref[1][2][0] = 1. + offset;
    mesh3D_ref[1][2][1] = 1.;
    mesh3D_ref[1][2][2] = 0.;
    mesh3D_ref[1][3][0] = 1.5 + offset;
    mesh3D_ref[1][3][1] = 0.;
    mesh3D_ref[1][3][2] = 1.;
    // Third cell
    mesh3D_ref[2].resize( 8 );
    mesh3D_ref[2][0][0] = 1.5 + offset;
    mesh3D_ref[2][0][1] = 0.;
    mesh3D_ref[2][0][2] = 0.;
    mesh3D_ref[2][1][0] = 3. + offset;
    mesh3D_ref[2][1][1] = 0.;
    mesh3D_ref[2][1][2] = 0.;
    mesh3D_ref[2][2][0] = 3. + offset;
    mesh3D_ref[2][2][1] = 1.;
    mesh3D_ref[2][2][2] = 0.;
    mesh3D_ref[2][3][0] = 2. + offset;
    mesh3D_ref[2][3][1] = 1.;
    mesh3D_ref[2][3][2] = 0.;
    mesh3D_ref[2][4][0] = 1.5 + offset;
    mesh3D_ref[2][4][1] = 0.;
    mesh3D_ref[2][4][2] = 1.;
    mesh3D_ref[2][5][0] = 3. + offset;
    mesh3D_ref[2][5][1] = 0.;
    mesh3D_ref[2][5][2] = 1.;
    mesh3D_ref[2][6][0] = 3. + offset;
    mesh3D_ref[2][6][1] = 1.;
    mesh3D_ref[2][6][2] = 1.;
    mesh3D_ref[2][7][0] = 1.5 + offset;
    mesh3D_ref[2][7][1] = 1.;
    mesh3D_ref[2][7][2] = 1.;
    // Fourth cell
    mesh3D_ref[3].resize( 8 );
    mesh3D_ref[3][0][0] = offset;
    mesh3D_ref[3][0][1] = 1.;
    mesh3D_ref[3][0][2] = 0.;
    mesh3D_ref[3][1][0] = 1. + offset;
    mesh3D_ref[3][1][1] = 1.;
    mesh3D_ref[3][1][2] = 0.;
    mesh3D_ref[3][2][0] = 1.5 + offset;
    mesh3D_ref[3][2][1] = 2.;
    mesh3D_ref[3][2][2] = 0.;
    mesh3D_ref[3][3][0] = offset;
    mesh3D_ref[3][3][1] = 2.;
    mesh3D_ref[3][3][2] = 0.;
    mesh3D_ref[3][4][0] = offset;
    mesh3D_ref[3][4][1] = 1.;
    mesh3D_ref[3][4][2] = 1.;
    mesh3D_ref[3][5][0] = 1.5 + offset;
    mesh3D_ref[3][5][1] = 1.;
    mesh3D_ref[3][5][2] = 1.;
    mesh3D_ref[3][6][0] = 1.5 + offset;
    mesh3D_ref[3][6][1] = 2.;
    mesh3D_ref[3][6][2] = 1.;
    mesh3D_ref[3][7][0] = offset;
    mesh3D_ref[3][7][1] = 2.;
    mesh3D_ref[3][7][2] = 1.;
    // Fifth cell
    mesh3D_ref[4].resize( 4 );
    mesh3D_ref[4][0][0] = 1. + offset;
    mesh3D_ref[4][0][1] = 1.;
    mesh3D_ref[4][0][2] = 0.;
    mesh3D_ref[4][1][0] = 2. + offset;
    mesh3D_ref[4][1][1] = 1.;
    mesh3D_ref[4][1][2] = 0.;
    mesh3D_ref[4][2][0] = 1.5 + offset;
    mesh3D_ref[4][2][1] = 2.;
    mesh3D_ref[4][2][2] = 0.;
    mesh3D_ref[4][3][0] = 1.5 + offset;
    mesh3D_ref[4][3][1] = 2.;
    mesh3D_ref[4][3][2] = 1.;
    // Sixth cell
    mesh3D_ref[5].resize( 8 );
    mesh3D_ref[5][0][0] = 2. + offset;
    mesh3D_ref[5][0][1] = 1.;
    mesh3D_ref[5][0][2] = 0.;
    mesh3D_ref[5][1][0] = 3. + offset;
    mesh3D_ref[5][1][1] = 1.;
    mesh3D_ref[5][1][2] = 0.;
    mesh3D_ref[5][2][0] = 3. + offset;
    mesh3D_ref[5][2][1] = 2.;
    mesh3D_ref[5][2][2] = 0.;
    mesh3D_ref[5][3][0] = 1.5 + offset;
    mesh3D_ref[5][3][1] = 2.;
    mesh3D_ref[5][3][2] = 0.;
    mesh3D_ref[5][4][0] = 1.5 + offset;
    mesh3D_ref[5][4][1] = 1.;
    mesh3D_ref[5][4][2] = 1.;
    mesh3D_ref[5][5][0] = 3. + offset;
    mesh3D_ref[5][5][1] = 1.;
    mesh3D_ref[5][5][2] = 1.;
    mesh3D_ref[5][6][0] = 3. + offset;
    mesh3D_ref[5][6][1] = 2.;
    mesh3D_ref[5][6][2] = 1.;
    mesh3D_ref[5][7][0] = 1.5 + offset;
    mesh3D_ref[5][7][1] = 2.;
    mesh3D_ref[5][7][2] = 1.;

    std::tie( cell_topologies_view, cells, coordinates ) =
        buildMixedMesh<DeviceType>( comm, dim_3 );
    TEST_EQUALITY( cell_topologies_view.extent( 0 ), n_cells );
    TEST_EQUALITY( cells.extent( 0 ), n_vertices );

    cell_topologies_view_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    Kokkos::deep_copy( cell_topologies_view_host, cell_topologies_view );
    cell_topology_ref = {
        {DTK_HEX_8, DTK_TET_4, DTK_HEX_8, DTK_HEX_8, DTK_TET_4, DTK_HEX_8}};
    for ( unsigned int i = 0; i < n_cells; ++i )
        TEST_EQUALITY( cell_topologies_view_host( i ), cell_topology_ref[i] );

    cells_host = Kokkos::create_mirror_view( cells );
    Kokkos::deep_copy( cells_host, cells );
    coordinates_host = Kokkos::create_mirror_view( coordinates );
    Kokkos::deep_copy( coordinates_host, coordinates );
    n = 0;
    for ( unsigned int i = 0; i < n_cells; ++i )
    {
        for ( unsigned int j = 0; j < mesh3D_ref[i].size(); ++j )
        {
            for ( unsigned int k = 0; k < dim_3; ++k )
            {
                unsigned int coord_pos = cells_host( n );
                TEST_EQUALITY( coordinates_host( coord_pos, k ),
                               mesh3D_ref[i][j][k] );
            }
            ++n;
        }
    }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MeshGenerator, simplex, DeviceType )
{
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();

    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<double **, DeviceType> coordinates;

    // 2D test
    std::vector<unsigned int> n_subdivisions = {{4, 3}};
    unsigned int n_cells = 2;
    for ( auto n_sub : n_subdivisions )
        n_cells *= n_sub;
    unsigned int constexpr dim_2 = 2;
    unsigned int constexpr n_vertices_per_tri = 3;
    int const comm_rank = comm->getRank();
    double offset = comm_rank * n_subdivisions[1];
    std::vector<std::array<std::array<double, dim_2>, n_vertices_per_tri>>
        tri_mesh_ref( n_cells );
    unsigned int n = 0;
    for ( unsigned int i = 0; i < n_subdivisions[1]; ++i )
        for ( unsigned int j = 0; j < n_subdivisions[0]; ++j )
        {
            tri_mesh_ref[n][0][0] = j;
            tri_mesh_ref[n][0][1] = i + offset;

            tri_mesh_ref[n][1][0] = j + 1;
            tri_mesh_ref[n][1][1] = i + offset;

            tri_mesh_ref[n][2][0] = j;
            tri_mesh_ref[n][2][1] = i + 1 + offset;

            ++n;

            tri_mesh_ref[n][0][0] = j + 1;
            tri_mesh_ref[n][0][1] = i + offset;

            tri_mesh_ref[n][1][0] = j + 1;
            tri_mesh_ref[n][1][1] = i + 1 + offset;

            tri_mesh_ref[n][2][0] = j;
            tri_mesh_ref[n][2][1] = i + 1 + offset;

            ++n;
        }

    std::tie( cell_topologies_view, cells, coordinates ) =
        buildSimplexMesh<DeviceType>( comm, n_subdivisions );
    TEST_EQUALITY( cell_topologies_view.extent( 0 ), n_cells );
    TEST_EQUALITY( cells.extent( 0 ), n_cells * n_vertices_per_tri );

    auto cell_topologies_view_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    Kokkos::deep_copy( cell_topologies_view_host, cell_topologies_view );
    for ( unsigned int i = 0; i < n_cells; ++i )
        TEST_EQUALITY( cell_topologies_view_host( i ), DTK_TRI_3 );

    auto cells_host = Kokkos::create_mirror_view( cells );
    Kokkos::deep_copy( cells_host, cells );
    auto coordinates_host = Kokkos::create_mirror_view( coordinates );
    Kokkos::deep_copy( coordinates_host, coordinates );
    n = 0;
    for ( unsigned int i = 0; i < n_cells; ++i )
    {
        for ( unsigned int j = 0; j < n_vertices_per_tri; ++j )
        {
            for ( unsigned int k = 0; k < dim_2; ++k )
            {
                unsigned int coord_pos = cells_host( n );
                TEST_EQUALITY( coordinates_host( coord_pos, k ),
                               tri_mesh_ref[i][j][k] );
            }

            ++n;
        }
    }

    // 3D test
    n_subdivisions = {{2, 2, 2}};
    n_cells = 5;
    for ( auto n_sub : n_subdivisions )
        n_cells *= n_sub;
    unsigned int constexpr n_vertices_per_tet = 4;
    offset = comm_rank * n_subdivisions[2];

    unsigned int constexpr dim_3 = 3;
    std::vector<std::array<std::array<double, dim_3>, n_vertices_per_tet>>
        tet_mesh_ref( n_cells );
    n = 0;
    for ( unsigned int i = 0; i < n_subdivisions[2]; i += 2 )
        for ( unsigned int j = 0; j < n_subdivisions[1]; ++j )
            for ( unsigned int k = 0; k < n_subdivisions[0]; ++k )
            {
                // First tet
                tet_mesh_ref[n][0][0] = k;
                tet_mesh_ref[n][0][1] = j;
                tet_mesh_ref[n][0][2] = i + offset;

                tet_mesh_ref[n][1][0] = k + 1;
                tet_mesh_ref[n][1][1] = j;
                tet_mesh_ref[n][1][2] = i + offset;

                tet_mesh_ref[n][2][0] = k;
                tet_mesh_ref[n][2][1] = j + 1;
                tet_mesh_ref[n][2][2] = i + offset;

                tet_mesh_ref[n][3][0] = k;
                tet_mesh_ref[n][3][1] = j;
                tet_mesh_ref[n][3][2] = i + 1 + offset;

                ++n;

                // Second tet
                tet_mesh_ref[n][0][0] = k + 1;
                tet_mesh_ref[n][0][1] = j;
                tet_mesh_ref[n][0][2] = i + offset;

                tet_mesh_ref[n][1][0] = k + 1;
                tet_mesh_ref[n][1][1] = j + 1;
                tet_mesh_ref[n][1][2] = i + offset;

                tet_mesh_ref[n][2][0] = k;
                tet_mesh_ref[n][2][1] = j + 1;
                tet_mesh_ref[n][2][2] = i + offset;

                tet_mesh_ref[n][3][0] = k + 1;
                tet_mesh_ref[n][3][1] = j + 1;
                tet_mesh_ref[n][3][2] = i + 1 + offset;

                ++n;

                // Third tet
                tet_mesh_ref[n][0][0] = k + 1;
                tet_mesh_ref[n][0][1] = j;
                tet_mesh_ref[n][0][2] = i + offset;

                tet_mesh_ref[n][1][0] = k + 1;
                tet_mesh_ref[n][1][1] = j + 1;
                tet_mesh_ref[n][1][2] = i + 1 + offset;

                tet_mesh_ref[n][2][0] = k;
                tet_mesh_ref[n][2][1] = j + 1;
                tet_mesh_ref[n][2][2] = i + offset;

                tet_mesh_ref[n][3][0] = k;
                tet_mesh_ref[n][3][1] = j;
                tet_mesh_ref[n][3][2] = i + 1 + offset;

                ++n;

                // Fourth tet
                tet_mesh_ref[n][0][0] = k;
                tet_mesh_ref[n][0][1] = j + 1;
                tet_mesh_ref[n][0][2] = i + offset;

                tet_mesh_ref[n][1][0] = k + 1;
                tet_mesh_ref[n][1][1] = j;
                tet_mesh_ref[n][1][2] = i + 1 + offset;

                tet_mesh_ref[n][2][0] = k + 1;
                tet_mesh_ref[n][2][1] = j + 1;
                tet_mesh_ref[n][2][2] = i + 1 + offset;

                tet_mesh_ref[n][3][0] = k;
                tet_mesh_ref[n][3][1] = j + 1;
                tet_mesh_ref[n][3][2] = i + 1 + offset;

                ++n;

                // Fifth tet
                tet_mesh_ref[n][0][0] = k + 1;
                tet_mesh_ref[n][0][1] = j;
                tet_mesh_ref[n][0][2] = i + offset;

                tet_mesh_ref[n][1][0] = k + 1;
                tet_mesh_ref[n][1][1] = j + 1;
                tet_mesh_ref[n][1][2] = i + 1 + offset;

                tet_mesh_ref[n][2][0] = k;
                tet_mesh_ref[n][2][1] = j;
                tet_mesh_ref[n][2][2] = i + 1 + offset;

                tet_mesh_ref[n][3][0] = k + 1;
                tet_mesh_ref[n][3][1] = j;
                tet_mesh_ref[n][3][2] = i + 1 + offset;

                ++n;

                // Sixth tet
                tet_mesh_ref[n][0][0] = k + 1;
                tet_mesh_ref[n][0][1] = j;
                tet_mesh_ref[n][0][2] = i + 1 + offset;

                tet_mesh_ref[n][1][0] = k + 1;
                tet_mesh_ref[n][1][1] = j + 1;
                tet_mesh_ref[n][1][2] = i + 1 + offset;

                tet_mesh_ref[n][2][0] = k;
                tet_mesh_ref[n][2][1] = j;
                tet_mesh_ref[n][2][2] = i + 1 + offset;

                tet_mesh_ref[n][3][0] = k + 1;
                tet_mesh_ref[n][3][1] = j;
                tet_mesh_ref[n][3][2] = i + 2 + offset;

                ++n;

                // Seventh tet
                tet_mesh_ref[n][0][0] = k + 1;
                tet_mesh_ref[n][0][1] = j + 1;
                tet_mesh_ref[n][0][2] = i + 1 + offset;

                tet_mesh_ref[n][1][0] = k;
                tet_mesh_ref[n][1][1] = j + 1;
                tet_mesh_ref[n][1][2] = i + 1 + offset;

                tet_mesh_ref[n][2][0] = k;
                tet_mesh_ref[n][2][1] = j;
                tet_mesh_ref[n][2][2] = i + 1 + offset;

                tet_mesh_ref[n][3][0] = k;
                tet_mesh_ref[n][3][1] = j + 1;
                tet_mesh_ref[n][3][2] = i + 2 + offset;

                ++n;

                // Eighth tet
                tet_mesh_ref[n][0][0] = k;
                tet_mesh_ref[n][0][1] = j;
                tet_mesh_ref[n][0][2] = i + 1 + offset;

                tet_mesh_ref[n][1][0] = k + 1;
                tet_mesh_ref[n][1][1] = j + 1;
                tet_mesh_ref[n][1][2] = i + 1 + offset;

                tet_mesh_ref[n][2][0] = k + 1;
                tet_mesh_ref[n][2][1] = j;
                tet_mesh_ref[n][2][2] = i + 2 + offset;

                tet_mesh_ref[n][3][0] = k;
                tet_mesh_ref[n][3][1] = j + 1;
                tet_mesh_ref[n][3][2] = i + 2 + offset;

                ++n;

                // Nineth tet
                tet_mesh_ref[n][0][0] = k;
                tet_mesh_ref[n][0][1] = j;
                tet_mesh_ref[n][0][2] = i + 2 + offset;

                tet_mesh_ref[n][1][0] = k;
                tet_mesh_ref[n][1][1] = j + 1;
                tet_mesh_ref[n][1][2] = i + 2 + offset;

                tet_mesh_ref[n][2][0] = k + 1;
                tet_mesh_ref[n][2][1] = j;
                tet_mesh_ref[n][2][2] = i + 2 + offset;

                tet_mesh_ref[n][3][0] = k;
                tet_mesh_ref[n][3][1] = j;
                tet_mesh_ref[n][3][2] = i + 1 + offset;

                ++n;

                // Tenth tet
                tet_mesh_ref[n][0][0] = k;
                tet_mesh_ref[n][0][1] = j + 1;
                tet_mesh_ref[n][0][2] = i + 2 + offset;

                tet_mesh_ref[n][1][0] = k + 1;
                tet_mesh_ref[n][1][1] = j + 1;
                tet_mesh_ref[n][1][2] = i + 2 + offset;

                tet_mesh_ref[n][2][0] = k + 1;
                tet_mesh_ref[n][2][1] = j;
                tet_mesh_ref[n][2][2] = i + 2 + offset;

                tet_mesh_ref[n][3][0] = k + 1;
                tet_mesh_ref[n][3][1] = j + 1;
                tet_mesh_ref[n][3][2] = i + 1 + offset;

                ++n;
            }

    std::tie( cell_topologies_view, cells, coordinates ) =
        buildSimplexMesh<DeviceType>( comm, n_subdivisions );
    TEST_EQUALITY( cell_topologies_view.extent( 0 ), n_cells );
    TEST_EQUALITY( cells.extent( 0 ), n_cells * n_vertices_per_tet );

    cell_topologies_view_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    Kokkos::deep_copy( cell_topologies_view_host, cell_topologies_view );
    for ( unsigned int i = 0; i < n_cells; ++i )
        TEST_EQUALITY( cell_topologies_view_host( i ), DTK_TET_4 );

    cells_host = Kokkos::create_mirror_view( cells );
    Kokkos::deep_copy( cells_host, cells );
    coordinates_host = Kokkos::create_mirror_view( coordinates );
    Kokkos::deep_copy( coordinates_host, coordinates );
    n = 0;
    for ( unsigned int i = 0; i < n_cells; ++i )
    {
        for ( unsigned int j = 0; j < tet_mesh_ref[i].size(); ++j )
        {
            for ( unsigned int k = 0; k < dim_3; ++k )
            {
                unsigned int coord_pos = cells_host( n );
                TEST_EQUALITY( coordinates_host( coord_pos, k ),
                               tet_mesh_ref[i][j][k] );
            }
            ++n;
        }
    }
}

// Include the test macros.
#include "DataTransferKitDiscretization_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    using DeviceType##NODE = typename NODE::device_type;                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MeshGenerator, structured,           \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MeshGenerator, mixed,                \
                                          DeviceType##NODE )                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MeshGenerator, simplex,              \
                                          DeviceType##NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
