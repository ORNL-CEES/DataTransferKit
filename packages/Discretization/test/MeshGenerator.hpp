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

#ifndef DTK_MESHGENERATOR_HPP
#define DTK_MESHGENERATOR_HPP

#include <DTK_CellTypes.h>
#include <DTK_DBC.hpp>
#include <DTK_Types.h>

#include <Kokkos_Core.hpp>

#include <mpi.h>

#include <vector>

// Compute the coordinates of the vertices of a simple 2D/3D slab domain
template <typename DeviceType>
Kokkos::View<DataTransferKit::Coordinate **, DeviceType>
computeCoordinates( unsigned int const n_vertices,
                    std::vector<unsigned int> const &n_subdivisions,
                    unsigned int comm_rank )
{
    unsigned int const dim = n_subdivisions.size();

    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> coordinates(
        "coordinates", n_vertices, dim );
    auto coordinates_host = Kokkos::create_mirror_view( coordinates );
    std::vector<unsigned int> current_vertex( dim, 0 );
    unsigned int proc_offset = n_subdivisions[dim - 1] * comm_rank;
    for ( unsigned int i = 0; i < n_vertices; ++i )
    {
        for ( unsigned int d = 0; d < dim; ++d )
        {
            unsigned int pos = current_vertex[d];
            unsigned int j = ( d == dim - 1 ) ? ( proc_offset + pos ) : pos;
            coordinates_host( i, d ) = static_cast<double>( j );
        }

        // Go to the next vertex
        ++current_vertex[0];
        if ( ( current_vertex[0] != 0 ) &&
             ( current_vertex[0] % ( n_subdivisions[0] + 1 ) ) == 0 )
        {
            current_vertex[0] = 0;
            ++current_vertex[1];
            if ( ( dim == 3 ) && ( current_vertex[1] != 0 ) &&
                 ( ( current_vertex[1] % ( n_subdivisions[1] + 1 ) ) == 0 ) )
            {
                current_vertex[1] = 0;
                ++current_vertex[2];
            }
        }
    }
    Kokkos::deep_copy( coordinates, coordinates_host );
    return coordinates;
}

template <typename DeviceType>
std::tuple<Kokkos::View<DTK_CellTopology *, DeviceType>,
           Kokkos::View<unsigned int *, DeviceType>,
           Kokkos::View<DataTransferKit::Coordinate **, DeviceType>>
buildStructuredMesh( MPI_Comm comm,
                     std::vector<unsigned int> const &n_subdivisions )
{
    // Build a 2D/3D structured mesh, i.e.,
    // ---------
    // | | | | |
    // ---------
    // | | | | |
    // ---------
    //
    // The mesh is offset on each rank in the latest direction (y in 2D and z in
    // 3D).

    unsigned int const dim = n_subdivisions.size();
    DTK_REQUIRE( ( dim == 2 ) || ( dim == 3 ) );

    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    unsigned int n_local_cells = 1;
    unsigned int n_vertices = 1;
    for ( auto n_sub : n_subdivisions )
    {
        n_local_cells *= n_sub;
        n_vertices *= n_sub + 1;
    }

    // Create the Kokkos::View of the coordinates
    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> coordinates =
        computeCoordinates<DeviceType>( n_vertices, n_subdivisions, comm_rank );

    // Create the Kokkos::View of the coordinates
    unsigned int const n_vertices_per_cell = std::pow( 2, dim );
    Kokkos::View<unsigned int *, DeviceType> cells(
        "cells", n_local_cells * n_vertices_per_cell );
    auto cells_host = Kokkos::create_mirror_view( cells );
    unsigned int n = 0;
    unsigned int const k_max = n_subdivisions[0];
    unsigned int const j_max = n_subdivisions[1];
    unsigned int const i_max = ( dim == 3 ) ? n_subdivisions[2] : 1;
    unsigned int const k_offset = k_max + 1;
    unsigned int const j_offset = j_max + 1;
    for ( unsigned int i = 0; i < i_max; ++i )
        for ( unsigned int j = 0; j < j_max; ++j )
            for ( unsigned int k = 0; k < k_max; ++k )
            {
                cells_host( n++ ) = k + j * k_offset + i * j_offset * k_offset;
                cells_host( n++ ) =
                    ( k + 1 ) + j * k_offset + i * j_offset * k_offset;
                cells_host( n++ ) =
                    ( k + 1 ) + ( j + 1 ) * k_offset + i * j_offset * k_offset;
                cells_host( n++ ) =
                    k + ( j + 1 ) * k_offset + i * j_offset * k_offset;

                if ( dim == 3 )
                {
                    cells_host( n++ ) =
                        k + j * k_offset + ( i + 1 ) * j_offset * k_offset;
                    cells_host( n++ ) = ( k + 1 ) + j * k_offset +
                                        ( i + 1 ) * j_offset * k_offset;
                    cells_host( n++ ) = ( k + 1 ) + ( j + 1 ) * k_offset +
                                        ( i + 1 ) * j_offset * k_offset;
                    cells_host( n++ ) = k + ( j + 1 ) * k_offset +
                                        ( i + 1 ) * j_offset * k_offset;
                }
            }
    Kokkos::deep_copy( cells, cells_host );

    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view(
        "cell_topologies", n_local_cells );
    Kokkos::deep_copy( cell_topologies_view,
                       ( dim == 2 ) ? DTK_QUAD_4 : DTK_HEX_8 );

    return std::make_tuple( cell_topologies_view, cells, coordinates );
}

template <typename DeviceType>
std::tuple<Kokkos::View<DTK_CellTopology *, DeviceType>,
           Kokkos::View<unsigned int *, DeviceType>,
           Kokkos::View<DataTransferKit::Coordinate **, DeviceType>>
buildMixedMesh( MPI_Comm comm, unsigned int const dim )
{
    // The mesh is composed of 4 quads + 2 triangles in 2D and 4 hex + 2 tets in
    // 3D. In 2D, the mesh looks like this:
    //
    // 7-----8-----9
    // |    / \    |
    // 3---4---5---6
    // |    \ /    |
    // 0-----1-----2
    //
    // The mesh is always offset in the x direction.

    DTK_REQUIRE( ( dim == 2 ) || ( dim == 3 ) );

    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    // Build the mesh
    unsigned int const n_local_hex_cells = 4;
    unsigned int const n_hex_vertices = ( dim == 2 ) ? 4 : 8;
    unsigned int const n_local_simplex_cells = 2;
    unsigned int const n_simplex_vertices = ( dim == 2 ) ? 3 : 4;
    unsigned int const n_local_cells =
        n_local_hex_cells + n_local_simplex_cells;
    unsigned int const offset_mesh = 3 * comm_rank;
    unsigned int const n_vertices = ( dim == 2 ) ? 10 : 19;

    // Create the Kokkos::View of the coordinates
    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> coordinates(
        "coordinates", n_vertices, dim );
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

    if ( dim == 3 )
    {
        for ( unsigned int i = 0; i < 10; ++i )
            coordinates_host( i, 2 ) = 0.;

        // Y=0 and Z=1
        coordinates_host( 10, 0 ) = offset_mesh;
        coordinates_host( 10, 1 ) = 0;
        coordinates_host( 10, 2 ) = 1;
        coordinates_host( 11, 0 ) = offset_mesh + 1.5;
        coordinates_host( 11, 1 ) = 0;
        coordinates_host( 11, 2 ) = 1;
        coordinates_host( 12, 0 ) = offset_mesh + 3;
        coordinates_host( 12, 1 ) = 0;
        coordinates_host( 12, 2 ) = 1;
        // Y=1 and Z=1
        coordinates_host( 13, 0 ) = offset_mesh;
        coordinates_host( 13, 1 ) = 1;
        coordinates_host( 13, 2 ) = 1;
        coordinates_host( 14, 0 ) = offset_mesh + 1.5;
        coordinates_host( 14, 1 ) = 1;
        coordinates_host( 14, 2 ) = 1;
        coordinates_host( 15, 0 ) = offset_mesh + 3;
        coordinates_host( 15, 1 ) = 1;
        coordinates_host( 15, 2 ) = 1;
        // Y=2 and Z=1
        coordinates_host( 16, 0 ) = offset_mesh;
        coordinates_host( 16, 1 ) = 2;
        coordinates_host( 16, 2 ) = 1;
        coordinates_host( 17, 0 ) = offset_mesh + 1.5;
        coordinates_host( 17, 1 ) = 2;
        coordinates_host( 17, 2 ) = 1;
        coordinates_host( 18, 0 ) = offset_mesh + 3;
        coordinates_host( 18, 1 ) = 2;
        coordinates_host( 18, 2 ) = 1;
    }

    Kokkos::deep_copy( coordinates, coordinates_host );

    // Create the Kokkos::View of the coordinates
    Kokkos::View<unsigned int *, DeviceType> cells(
        "cells", n_hex_vertices * n_local_hex_cells +
                     n_simplex_vertices * n_local_simplex_cells );
    auto cells_host = Kokkos::create_mirror_view( cells );
    unsigned int n = 0;
    // First cell
    /*
     *  --
     *  | \
     *  ---
     */
    cells_host( n++ ) = 0;
    cells_host( n++ ) = 1;
    cells_host( n++ ) = 4;
    cells_host( n++ ) = 3;
    if ( dim == 3 )
    {
        cells_host( n++ ) = 10;
        cells_host( n++ ) = 11;
        cells_host( n++ ) = 14;
        cells_host( n++ ) = 13;
    }
    // Second cell
    // --
    // \/
    cells_host( n++ ) = 1;
    cells_host( n++ ) = 5;
    cells_host( n++ ) = 4;
    if ( dim == 3 )
    {
        cells_host( n++ ) = 11;
    }
    // Third cell
    //  --
    // / |
    // ---
    cells_host( n++ ) = 1;
    cells_host( n++ ) = 2;
    cells_host( n++ ) = 6;
    cells_host( n++ ) = 5;
    if ( dim == 3 )
    {
        cells_host( n++ ) = 11;
        cells_host( n++ ) = 12;
        cells_host( n++ ) = 15;
        cells_host( n++ ) = 14;
    }
    // Fourth cell
    // ---
    // | /
    // --
    cells_host( n++ ) = 3;
    cells_host( n++ ) = 4;
    cells_host( n++ ) = 8;
    cells_host( n++ ) = 7;
    if ( dim == 3 )
    {
        cells_host( n++ ) = 13;
        cells_host( n++ ) = 14;
        cells_host( n++ ) = 17;
        cells_host( n++ ) = 16;
    }
    // Fifth cell
    /*
     * /\
     * --
     */
    cells_host( n++ ) = 4;
    cells_host( n++ ) = 5;
    cells_host( n++ ) = 8;
    if ( dim == 3 )
    {
        cells_host( n++ ) = 17;
    }
    // Sixth cell
    // ---
    // \ |
    //  --
    cells_host( n++ ) = 5;
    cells_host( n++ ) = 6;
    cells_host( n++ ) = 9;
    cells_host( n++ ) = 8;
    if ( dim == 3 )
    {
        cells_host( n++ ) = 14;
        cells_host( n++ ) = 15;
        cells_host( n++ ) = 18;
        cells_host( n++ ) = 17;
    }
    Kokkos::deep_copy( cells, cells_host );

    // Create the Kokkos::View of the topologies
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view(
        "cell_topologies", n_local_cells );
    auto cell_topologies_view_host =
        Kokkos::create_mirror_view( cell_topologies_view );
    cell_topologies_view_host( 0 ) = ( dim == 2 ) ? DTK_QUAD_4 : DTK_HEX_8;
    cell_topologies_view_host( 1 ) = ( dim == 2 ) ? DTK_TRI_3 : DTK_TET_4;
    cell_topologies_view_host( 2 ) = ( dim == 2 ) ? DTK_QUAD_4 : DTK_HEX_8;
    cell_topologies_view_host( 3 ) = ( dim == 2 ) ? DTK_QUAD_4 : DTK_HEX_8;
    cell_topologies_view_host( 4 ) = ( dim == 2 ) ? DTK_TRI_3 : DTK_TET_4;
    cell_topologies_view_host( 5 ) = ( dim == 2 ) ? DTK_QUAD_4 : DTK_HEX_8;
    Kokkos::deep_copy( cell_topologies_view, cell_topologies_view_host );

    return std::make_tuple( cell_topologies_view, cells, coordinates );
}

template <typename DeviceType>
std::tuple<Kokkos::View<DTK_CellTopology *, DeviceType>,
           Kokkos::View<unsigned int *, DeviceType>,
           Kokkos::View<DataTransferKit::Coordinate **, DeviceType>>
buildSimplexMesh( MPI_Comm comm, std::vector<unsigned int> &n_subdivisions )
{
    // The mesh looks like this in 2D
    // -----------
    // |\|\|\|\|\|
    // -----------
    // |\|\|\|\|\|
    // -----------
    //
    // In 3D, the mesh is composed of hexahedra that have been divided in five
    // tetrahedra see
    // http://www.matematicasvisuales.com/english/html/geometry/space/voltetra.html
    //
    // The mesh is offset on each rank in the latest direction (y in 2D and z in
    // 3D).

    unsigned int const dim = n_subdivisions.size();
    DTK_REQUIRE( ( dim == 2 ) || ( dim == 3 ) );

    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );
    unsigned int n_local_cells = ( dim == 2 ) ? 2 : 5;
    unsigned int n_vertices = 1;
    if ( ( dim == 3 ) && ( n_subdivisions.back() % 2 == 1 ) )
        ++n_subdivisions.back();
    for ( auto n_sub : n_subdivisions )
    {
        n_local_cells *= n_sub;
        n_vertices *= n_sub + 1;
    }

    // Create the Kokkos::View of the coordinates
    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> coordinates =
        computeCoordinates<DeviceType>( n_vertices, n_subdivisions, comm_rank );

    unsigned int const n_vertices_per_cell = ( dim == 2 ) ? 3 : 4;
    Kokkos::View<unsigned int *, DeviceType> cells(
        "cells", n_local_cells * n_vertices_per_cell );
    auto cells_host = Kokkos::create_mirror_view( cells );
    unsigned int n = 0;
    unsigned int const k_max = n_subdivisions[0];
    unsigned int const j_max = n_subdivisions[1];
    unsigned int const i_max = ( dim == 2 ) ? 0 : n_subdivisions[2];
    unsigned int const k_offset = k_max + 1;
    unsigned int const j_offset = j_max + 1;
    if ( dim == 2 )
    {
        for ( unsigned int j = 0; j < j_max; ++j )
            for ( unsigned int k = 0; k < k_max; ++k )
            {
                // First simplex
                /*
                 *  |\
                 *  --
                 */
                cells_host( n++ ) = k + j * k_offset;
                cells_host( n++ ) = ( k + 1 ) + j * k_offset;
                cells_host( n++ ) = k + ( j + 1 ) * k_offset;

                // Second simplex
                // --
                // \|
                cells_host( n++ ) = ( k + 1 ) + j * k_offset;
                cells_host( n++ ) = ( k + 1 ) + ( j + 1 ) * k_offset;
                cells_host( n++ ) = k + ( j + 1 ) * k_offset;
            }
    }
    else
    {
        // In 3D, to ensure the mesh is confirming, we will divided each
        // subdivision cube into two hexahedra composed of five tetrahedra each
        for ( unsigned int i = 0; i < i_max; i += 2 )
            for ( unsigned int j = 0; j < j_max; ++j )
                for ( unsigned int k = 0; k < k_max; ++k )
                {
                    // First hexahedra
                    // First simplex
                    cells_host( n++ ) =
                        k + j * k_offset + i * k_offset * j_offset;
                    cells_host( n++ ) =
                        ( k + 1 ) + j * k_offset + i * k_offset * j_offset;
                    cells_host( n++ ) =
                        k + ( j + 1 ) * k_offset + i * k_offset * j_offset;
                    cells_host( n++ ) =
                        k + j * k_offset + ( i + 1 ) * k_offset * j_offset;
                    // Second simplex
                    cells_host( n++ ) =
                        ( k + 1 ) + j * k_offset + i * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + ( j + 1 ) * k_offset +
                                        i * k_offset * j_offset;
                    cells_host( n++ ) =
                        k + ( j + 1 ) * k_offset + i * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + ( j + 1 ) * k_offset +
                                        ( i + 1 ) * k_offset * j_offset;
                    // Third simplex
                    cells_host( n++ ) =
                        ( k + 1 ) + j * k_offset + i * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + ( j + 1 ) * k_offset +
                                        ( i + 1 ) * k_offset * j_offset;
                    cells_host( n++ ) =
                        k + ( j + 1 ) * k_offset + i * k_offset * j_offset;
                    cells_host( n++ ) =
                        k + j * k_offset + ( i + 1 ) * k_offset * j_offset;
                    // Fourth simplex
                    cells_host( n++ ) =
                        k + ( j + 1 ) * k_offset + i * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + j * k_offset +
                                        ( i + 1 ) * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + ( j + 1 ) * k_offset +
                                        ( i + 1 ) * k_offset * j_offset;
                    cells_host( n++ ) = k + ( j + 1 ) * k_offset +
                                        ( i + 1 ) * k_offset * j_offset;
                    // Fifth simplex
                    cells_host( n++ ) =
                        ( k + 1 ) + j * k_offset + i * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + ( j + 1 ) * k_offset +
                                        ( i + 1 ) * k_offset * j_offset;
                    cells_host( n++ ) =
                        k + j * k_offset + ( i + 1 ) * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + j * k_offset +
                                        ( i + 1 ) * k_offset * j_offset;

                    // Second hexahedra on the top of the first one. The
                    // hexahedra is rotated of pi/2 to make the mesh conforming.
                    // First simplex
                    cells_host( n++ ) = ( k + 1 ) + j * k_offset +
                                        ( i + 1 ) * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + ( j + 1 ) * k_offset +
                                        ( i + 1 ) * k_offset * j_offset;
                    cells_host( n++ ) =
                        k + j * k_offset + ( i + 1 ) * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + j * k_offset +
                                        ( i + 2 ) * k_offset * j_offset;
                    // Second simplex
                    cells_host( n++ ) = ( k + 1 ) + ( j + 1 ) * k_offset +
                                        ( i + 1 ) * k_offset * j_offset;
                    cells_host( n++ ) = k + ( j + 1 ) * k_offset +
                                        ( i + 1 ) * k_offset * j_offset;
                    cells_host( n++ ) =
                        k + j * k_offset + ( i + 1 ) * k_offset * j_offset;
                    cells_host( n++ ) = k + ( j + 1 ) * k_offset +
                                        ( i + 2 ) * k_offset * j_offset;
                    // Third simplex
                    cells_host( n++ ) =
                        k + j * k_offset + ( i + 1 ) * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + ( j + 1 ) * k_offset +
                                        ( i + 1 ) * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + j * k_offset +
                                        ( i + 2 ) * k_offset * j_offset;
                    cells_host( n++ ) = k + ( j + 1 ) * k_offset +
                                        ( i + 2 ) * k_offset * j_offset;
                    // Fourth simplex
                    cells_host( n++ ) =
                        k + j * k_offset + ( i + 2 ) * k_offset * j_offset;
                    cells_host( n++ ) = k + ( j + 1 ) * k_offset +
                                        ( i + 2 ) * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + j * k_offset +
                                        ( i + 2 ) * k_offset * j_offset;
                    cells_host( n++ ) =
                        k + j * k_offset + ( i + 1 ) * k_offset * j_offset;
                    // Fifth simplex
                    cells_host( n++ ) = k + ( j + 1 ) * k_offset +
                                        ( i + 2 ) * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + ( j + 1 ) * k_offset +
                                        ( i + 2 ) * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + j * k_offset +
                                        ( i + 2 ) * k_offset * j_offset;
                    cells_host( n++ ) = ( k + 1 ) + ( j + 1 ) * k_offset +
                                        ( i + 1 ) * k_offset * j_offset;
                }
    }
    Kokkos::deep_copy( cells, cells_host );

    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view(
        "cell_topologies", n_local_cells );
    Kokkos::deep_copy( cell_topologies_view,
                       ( dim == 2 ) ? DTK_TRI_3 : DTK_TET_4 );

    return std::make_tuple( cell_topologies_view, cells, coordinates );
}

#endif
