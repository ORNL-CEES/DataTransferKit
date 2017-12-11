/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#include "MeshGenerator.hpp"
#include <DTK_PointSearch.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_UnitTestHarness.hpp>

template <typename DeviceType>
Kokkos::View<double * [3], DeviceType>
getPointsCoord3D( Teuchos::RCP<const Teuchos::Comm<int>> comm )
{
    int const comm_rank = comm->getRank();
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
    else if ( comm_rank == 1 )
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

// The `out` and `success` parameters come from the Teuchos unit testing macros
// expansion.
template <int dim, typename DeviceType>
void checkReferencePoints(
    Kokkos::View<int *, DeviceType> ranks,
    Kokkos::View<int *, DeviceType> cell_indices,
    Kokkos::View<DataTransferKit::Point *, DeviceType> reference_points,
    Kokkos::View<unsigned int *, DeviceType> offset_view,
    std::vector<
        std::vector<std::tuple<int, int, std::array<double, dim>>>> const
        &ref_sol,
    bool &success, Teuchos::FancyOStream &out )
{
    auto ranks_host = Kokkos::create_mirror_view( ranks );
    Kokkos::deep_copy( ranks_host, ranks );
    auto cell_indices_host = Kokkos::create_mirror_view( cell_indices );
    Kokkos::deep_copy( cell_indices_host, cell_indices );
    auto reference_points_host = Kokkos::create_mirror_view( reference_points );
    Kokkos::deep_copy( reference_points_host, reference_points );
    auto offset_view_host = Kokkos::create_mirror_view( offset_view );
    Kokkos::deep_copy( offset_view_host, offset_view );
    for ( unsigned int i = 0; i < offset_view_host.size() - 1; ++i )
        for ( unsigned int j = offset_view_host( i );
              j < offset_view_host( i + 1 ); ++j )
        {
            int rank = ranks_host( j );
            int cell_index = cell_indices_host( j );
            DataTransferKit::Point point = reference_points_host( j );
            bool pt_found = false;
            for ( unsigned int k = 0; k < ref_sol[i].size(); ++k )
            {
                auto const &ref_query = ref_sol[i][k];
                if ( ( rank == std::get<0>( ref_query ) ) &&
                     ( cell_index == std::get<1>( ref_query ) ) )
                {
                    bool same_coord = true;
                    for ( unsigned int d = 0; d < dim; ++d )
                        if ( std::abs( std::get<2>( ref_query )[d] -
                                       point[d] ) > 1e-14 )
                            same_coord = false;
                    if ( same_coord )
                    {
                        pt_found = true;
                        break;
                    }
                }
            }
            TEST_EQUALITY( pt_found, true );
        }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( PointSearch, one_topo_three_dim, DeviceType )
{
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int const comm_rank = comm->getRank();
    unsigned int constexpr dim = 3;
    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view;
    Kokkos::View<unsigned int *, DeviceType> cells;
    Kokkos::View<double **, DeviceType> coordinates;
    Kokkos::View<double * [3], DeviceType> points_coord;
    std::vector<unsigned int> n_subdivisions = {{5, 5, 3}};
    std::tie( cell_topologies_view, cells, coordinates ) =
        buildStructuredMesh<DeviceType>( comm, n_subdivisions );
    points_coord = getPointsCoord3D<DeviceType>( comm );

    // We are now done with building the mesh and we can do the search

    DataTransferKit::PointSearch<DeviceType> pt_search(
        comm, cell_topologies_view, cells, coordinates, points_coord, "all" );

    Kokkos::View<int *, DeviceType> ranks;
    Kokkos::View<int *, DeviceType> cell_indices;
    Kokkos::View<DataTransferKit::Point *, DeviceType> reference_points;
    Kokkos::View<unsigned int *, DeviceType> offset_view;
    std::tie( ranks, cell_indices, reference_points, offset_view ) =
        pt_search.getSearchResults();

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
    std::vector<std::vector<std::tuple<int, int, PtCoord>>> ref_sol;
    if ( comm_rank == 0 )
    {
        ref_sol.resize( 5 );

        // First query
        PtCoord ref_frame_1 = {{0., 0., 0.}};
        std::vector<std::tuple<int, int, PtCoord>> queries_1( 1 );
        queries_1[0] = std::make_tuple( 0, 0, ref_frame_1 );
        ref_sol[0] = queries_1;

        // Second query
        PtCoord ref_frame_2 = {{-0.5, 0.5, -0.5}};
        std::vector<std::tuple<int, int, PtCoord>> queries_2( 1 );
        queries_2[0] = std::make_tuple( 1, 11, ref_frame_2 );
        ref_sol[1] = queries_2;

        // Third query
        PtCoord ref_frame_3_0 = {{0.5, -1, -0.5}};
        PtCoord ref_frame_3_1 = {{0.5, 1, -0.5}};
        std::vector<std::tuple<int, int, PtCoord>> queries_3( 2 );
        queries_3[0] = std::make_tuple( 1, 12, ref_frame_3_0 );
        queries_3[1] = std::make_tuple( 1, 7, ref_frame_3_1 );
        ref_sol[2] = queries_3;

        // Fourth query
        PtCoord ref_frame_4_0 = {{0., -1., -1.}};
        PtCoord ref_frame_4_1 = {{0., 1., -1.}};
        PtCoord ref_frame_4_2 = {{0., 1., 1.}};
        PtCoord ref_frame_4_3 = {{0., -1., 1}};
        std::vector<std::tuple<int, int, PtCoord>> queries_4( 4 );
        queries_4[0] = std::make_tuple( 1, 12, ref_frame_4_0 );
        queries_4[1] = std::make_tuple( 1, 7, ref_frame_4_1 );
        queries_4[2] = std::make_tuple( 0, 57, ref_frame_4_2 );
        queries_4[3] = std::make_tuple( 0, 62, ref_frame_4_3 );
        ref_sol[3] = queries_4;

        // Fifth query
        PtCoord ref_frame_5_0 = {{-1., -1., 1.}};
        PtCoord ref_frame_5_1 = {{-1., 1., 1}};
        PtCoord ref_frame_5_2 = {{1., -1., 1.}};
        PtCoord ref_frame_5_3 = {{1., 1., 1.}};
        PtCoord ref_frame_5_4 = {{-1., -1., -1.}};
        PtCoord ref_frame_5_5 = {{-1., 1, -1.}};
        PtCoord ref_frame_5_6 = {{1, -1., -1.}};
        PtCoord ref_frame_5_7 = {{1., 1., -1.}};
        std::vector<std::tuple<int, int, PtCoord>> queries_5( 8 );
        queries_5[0] = std::make_tuple( 0, 37, ref_frame_5_0 );
        queries_5[1] = std::make_tuple( 0, 32, ref_frame_5_1 );
        queries_5[2] = std::make_tuple( 0, 36, ref_frame_5_2 );
        queries_5[3] = std::make_tuple( 0, 31, ref_frame_5_3 );
        queries_5[4] = std::make_tuple( 0, 62, ref_frame_5_4 );
        queries_5[5] = std::make_tuple( 0, 57, ref_frame_5_5 );
        queries_5[6] = std::make_tuple( 0, 61, ref_frame_5_6 );
        queries_5[7] = std::make_tuple( 0, 56, ref_frame_5_7 );
        ref_sol[4] = queries_5;
    }
    else if ( comm_rank == 1 )
    {
        ref_sol.resize( 5 );

        // First query
        PtCoord ref_frame_1 = {{0., 0., 0.}};
        std::vector<std::tuple<int, int, PtCoord>> queries_1( 1 );
        queries_1[0] = std::make_tuple( 0, 31, ref_frame_1 );
        ref_sol[0] = queries_1;

        // Second query
        PtCoord ref_frame_2 = {{-0.5, 0.5, -0.5}};
        std::vector<std::tuple<int, int, PtCoord>> queries_2( 1 );
        queries_2[0] = std::make_tuple( 1, 42, ref_frame_2 );
        ref_sol[1] = queries_2;

        // Third query
        PtCoord ref_frame_3_0 = {{0.5, -1, -0.5}};
        PtCoord ref_frame_3_1 = {{0.5, 1, -0.5}};
        std::vector<std::tuple<int, int, PtCoord>> queries_3( 2 );
        queries_3[0] = std::make_tuple( 1, 43, ref_frame_3_0 );
        queries_3[1] = std::make_tuple( 1, 38, ref_frame_3_1 );
        ref_sol[2] = queries_3;

        // Fourth query
        PtCoord ref_frame_4_0 = {{0., -1., -1.}};
        PtCoord ref_frame_4_1 = {{0., 1., -1.}};
        PtCoord ref_frame_4_2 = {{0., 1., 1.}};
        PtCoord ref_frame_4_3 = {{0., -1., 1}};
        std::vector<std::tuple<int, int, PtCoord>> queries_4( 4 );
        queries_4[0] = std::make_tuple( 1, 43, ref_frame_4_0 );
        queries_4[1] = std::make_tuple( 1, 38, ref_frame_4_1 );
        queries_4[2] = std::make_tuple( 1, 13, ref_frame_4_2 );
        queries_4[3] = std::make_tuple( 1, 18, ref_frame_4_3 );
        ref_sol[3] = queries_4;

        // Fifth query
        PtCoord ref_frame_5_0 = {{-1., -1., 1.}};
        PtCoord ref_frame_5_1 = {{-1., 1., 1}};
        PtCoord ref_frame_5_2 = {{1., -1., 1.}};
        PtCoord ref_frame_5_3 = {{1., 1., 1.}};
        PtCoord ref_frame_5_4 = {{-1., -1., -1.}};
        PtCoord ref_frame_5_5 = {{-1., 1, -1.}};
        PtCoord ref_frame_5_6 = {{1, -1., -1.}};
        PtCoord ref_frame_5_7 = {{1., 1., -1.}};
        std::vector<std::tuple<int, int, PtCoord>> queries_5( 8 );
        queries_5[0] = std::make_tuple( 0, 68, ref_frame_5_0 );
        queries_5[1] = std::make_tuple( 0, 63, ref_frame_5_1 );
        queries_5[2] = std::make_tuple( 0, 67, ref_frame_5_2 );
        queries_5[3] = std::make_tuple( 0, 62, ref_frame_5_3 );
        queries_5[4] = std::make_tuple( 1, 18, ref_frame_5_4 );
        queries_5[5] = std::make_tuple( 1, 13, ref_frame_5_5 );
        queries_5[6] = std::make_tuple( 1, 17, ref_frame_5_6 );
        queries_5[7] = std::make_tuple( 1, 12, ref_frame_5_7 );
        ref_sol[4] = queries_5;
    }

    // Check the results
    checkReferencePoints<dim, DeviceType>( ranks, cell_indices,
                                           reference_points, offset_view,
                                           ref_sol, success, out );
}

// Include the test macros.
#include "DataTransferKitDiscretization_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    using DeviceType##NODE = typename NODE::device_type;                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( PointSearch, one_topo_three_dim,     \
                                          DeviceType##NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
