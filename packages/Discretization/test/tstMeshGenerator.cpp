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
#include <Teuchos_UnitTestHarness.hpp>

#include <array>
#include <fstream>

std::tuple<std::vector<std::vector<DataTransferKit::Coordinate>>,
           std::vector<unsigned int>>
readInputFile( std::string const &filename )
{
    std::ifstream file( filename );
    DTK_REQUIRE( file.is_open() );

    unsigned int dim = 0;
    file >> dim;

    // Read the coordinates of the vertices
    unsigned int n_vertices = 0;
    file >> n_vertices;
    std::vector<std::vector<DataTransferKit::Coordinate>> coordinates_ref(
        n_vertices, std::vector<DataTransferKit::Coordinate>( dim, 0. ) );
    for ( unsigned int i = 0; i < n_vertices; ++i )
        for ( unsigned int j = 0; j < dim; ++j )
        {
            file >> coordinates_ref[i][j];
        }

    // Read the vertex IDs associated to each cell.
    unsigned int n_cells = 0;
    file >> n_cells;
    // We do know not the size of cell_ref because different cells can have a
    // different number of vertices
    std::vector<unsigned int> cells_ref;
    for ( unsigned int i = 0; i < n_cells; ++i )
    {
        unsigned int n_vertex_per_cell = 0;
        file >> n_vertex_per_cell;
        for ( unsigned int j = 0; j < n_vertex_per_cell; ++j )
        {
            unsigned int val = 0;
            file >> val;
            cells_ref.push_back( val );
        }
    }

    file.close();

    return std::make_tuple( coordinates_ref, cells_ref );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MeshGenerator, structured, DeviceType )
{
    MPI_Comm comm = MPI_COMM_WORLD;
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> coordinates;

    // 2D test
    std::string filename = "structured_2d.txt";
    std::vector<std::vector<DataTransferKit::Coordinate>> coordinates_ref;
    std::vector<unsigned int> cells_ref;
    std::tie( coordinates_ref, cells_ref ) = readInputFile( filename );
    // Move mesh according to the rank
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    std::vector<unsigned int> n_subdivisions = {{4, 3}};
    double offset = n_subdivisions[1] * comm_rank;
    for ( auto &coord : coordinates_ref )
        coord[1] += offset;

    std::tie( cell_topologies_view, cells, coordinates ) =
        buildStructuredMesh<DeviceType>( comm, n_subdivisions );

    // Check view size
    unsigned int n_vertices = coordinates_ref.size();
    unsigned int n_cells = 1;
    for ( auto n_sub : n_subdivisions )
        n_cells *= n_sub;
    TEST_EQUALITY( cell_topologies_view.extent( 0 ), n_cells );
    TEST_EQUALITY( cells.extent( 0 ), cells_ref.size() );
    TEST_EQUALITY( coordinates.extent( 0 ), n_vertices );

    // Check topology
    auto cell_topologies_view_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    Kokkos::deep_copy( cell_topologies_view_host, cell_topologies_view );
    for ( unsigned int i = 0; i < n_cells; ++i )
        TEST_EQUALITY( cell_topologies_view_host( i ), DTK_QUAD_4 );

    // Check cells
    auto cells_host = Kokkos::create_mirror_view( cells );
    Kokkos::deep_copy( cells_host, cells );
    TEST_COMPARE_ARRAYS( cells_host, cells_ref );

    // Check coordinates
    unsigned int dim = 2;
    auto coordinates_host = Kokkos::create_mirror_view( coordinates );
    Kokkos::deep_copy( coordinates_host, coordinates );
    for ( unsigned int i = 0; i < n_vertices; ++i )
        for ( unsigned int j = 0; j < dim; ++j )
            TEST_EQUALITY( coordinates_host( i, j ), coordinates_ref[i][j] );

    // 3D test
    filename = "structured_3d.txt";
    std::tie( coordinates_ref, cells_ref ) = readInputFile( filename );
    // Move mesh according to the rank
    n_subdivisions = {{2, 3, 4}};
    offset = n_subdivisions[2] * comm_rank;
    for ( auto &coord : coordinates_ref )
        coord[2] += offset;

    std::tie( cell_topologies_view, cells, coordinates ) =
        buildStructuredMesh<DeviceType>( comm, n_subdivisions );

    n_vertices = coordinates_ref.size();
    n_cells = 1;
    for ( auto n_sub : n_subdivisions )
        n_cells *= n_sub;

    TEST_EQUALITY( cell_topologies_view.extent( 0 ), n_cells );
    TEST_EQUALITY( cells.extent( 0 ), cells_ref.size() );
    TEST_EQUALITY( coordinates.extent( 0 ), n_vertices );

    // Check topology
    cell_topologies_view_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    Kokkos::deep_copy( cell_topologies_view_host, cell_topologies_view );
    for ( unsigned int i = 0; i < n_cells; ++i )
        TEST_EQUALITY( cell_topologies_view_host( i ), DTK_HEX_8 );

    // Check cells
    cells_host = Kokkos::create_mirror_view( cells );
    Kokkos::deep_copy( cells_host, cells );
    TEST_COMPARE_ARRAYS( cells_host, cells_ref );

    // Check coordinates
    dim = 3;
    coordinates_host = Kokkos::create_mirror_view( coordinates );
    Kokkos::deep_copy( coordinates_host, coordinates );
    for ( unsigned int i = 0; i < n_vertices; ++i )
        for ( unsigned int j = 0; j < dim; ++j )
            TEST_EQUALITY( coordinates_host( i, j ), coordinates_ref[i][j] );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MeshGenerator, mixed, DeviceType )
{
    MPI_Comm comm = MPI_COMM_WORLD;

    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> coordinates;

    // 2D test
    std::string filename = "mixed_2d.txt";
    std::vector<std::vector<DataTransferKit::Coordinate>> coordinates_ref;
    std::vector<unsigned int> cells_ref;
    std::tie( coordinates_ref, cells_ref ) = readInputFile( filename );
    unsigned int dim = 2;
    // Move mesh according to the rank
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    double offset = 3. * comm_rank;
    for ( auto &coord : coordinates_ref )
        coord[0] += offset;

    std::tie( cell_topologies_view, cells, coordinates ) =
        buildMixedMesh<DeviceType>( comm, dim );

    // Check view size
    unsigned int n_vertices = coordinates_ref.size();
    unsigned int n_cells = 6;
    TEST_EQUALITY( cell_topologies_view.extent( 0 ), n_cells );
    TEST_EQUALITY( cells.extent( 0 ), cells_ref.size() );
    TEST_EQUALITY( coordinates.extent( 0 ), n_vertices );

    // Check topology
    auto cell_topologies_view_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    Kokkos::deep_copy( cell_topologies_view_host, cell_topologies_view );
    std::vector<DTK_CellTopology> cell_topology_ref = {
        {DTK_QUAD_4, DTK_TRI_3, DTK_QUAD_4, DTK_QUAD_4, DTK_TRI_3, DTK_QUAD_4}};
    for ( unsigned int i = 0; i < n_cells; ++i )
        TEST_EQUALITY( cell_topologies_view_host( i ), cell_topology_ref[i] );

    // Check cells
    auto cells_host = Kokkos::create_mirror_view( cells );
    Kokkos::deep_copy( cells_host, cells );
    TEST_COMPARE_ARRAYS( cells_host, cells_ref );

    // Check coordinates
    auto coordinates_host = Kokkos::create_mirror_view( coordinates );
    Kokkos::deep_copy( coordinates_host, coordinates );
    for ( unsigned int i = 0; i < n_vertices; ++i )
        for ( unsigned int j = 0; j < dim; ++j )
            TEST_EQUALITY( coordinates_host( i, j ), coordinates_ref[i][j] );

    // 3D test
    filename = "mixed_3d.txt";
    std::tie( coordinates_ref, cells_ref ) = readInputFile( filename );
    dim = 3;
    // Move mesh according to the rank
    for ( auto &coord : coordinates_ref )
        coord[0] += offset;

    std::tie( cell_topologies_view, cells, coordinates ) =
        buildMixedMesh<DeviceType>( comm, dim );

    // Check view size
    n_vertices = coordinates_ref.size();
    TEST_EQUALITY( cell_topologies_view.extent( 0 ), n_cells );
    TEST_EQUALITY( cells.extent( 0 ), cells_ref.size() );
    TEST_EQUALITY( coordinates.extent( 0 ), n_vertices );

    // Check topology
    Kokkos::deep_copy( cell_topologies_view_host, cell_topologies_view );
    cell_topology_ref = {
        {DTK_HEX_8, DTK_TET_4, DTK_HEX_8, DTK_HEX_8, DTK_TET_4, DTK_HEX_8}};
    for ( unsigned int i = 0; i < n_cells; ++i )
        TEST_EQUALITY( cell_topologies_view_host( i ), cell_topology_ref[i] );

    // Check cells
    cells_host = Kokkos::create_mirror_view( cells );
    Kokkos::deep_copy( cells_host, cells );
    TEST_COMPARE_ARRAYS( cells_host, cells_ref );

    // Check coordinates
    coordinates_host = Kokkos::create_mirror_view( coordinates );
    Kokkos::deep_copy( coordinates_host, coordinates );
    for ( unsigned int i = 0; i < n_vertices; ++i )
        for ( unsigned int j = 0; j < dim; ++j )
            TEST_EQUALITY( coordinates_host( i, j ), coordinates_ref[i][j] );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MeshGenerator, simplex, DeviceType )
{
    MPI_Comm comm = MPI_COMM_WORLD;

    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> coordinates;

    // 2D test
    std::vector<unsigned int> n_subdivisions = {{4, 3}};
    unsigned int n_cells = 2;
    for ( auto n_sub : n_subdivisions )
        n_cells *= n_sub;
    unsigned int constexpr dim_2 = 2;
    unsigned int constexpr n_vertices_per_tri = 3;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    double offset = comm_rank * n_subdivisions[1];
    std::vector<std::array<std::array<DataTransferKit::Coordinate, dim_2>,
                           n_vertices_per_tri>>
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
    std::vector<std::array<std::array<DataTransferKit::Coordinate, dim_3>,
                           n_vertices_per_tet>>
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
#include "DataTransferKit_ETIHelperMacros.h"

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
