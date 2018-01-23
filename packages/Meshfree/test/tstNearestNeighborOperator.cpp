/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#include <Teuchos_UnitTestHarness.hpp>

#include <DTK_DBC.hpp>                              // DataTransferKitException
#include <DTK_DetailsDistributedSearchTreeImpl.hpp> // epsilon
#include <DTK_NearestNeighborOperator.hpp>
#include <Kokkos_Core.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Distributor.hpp>
#include <Tpetra_Map.hpp>

#include <array>
#include <numeric>
#include <random>
#include <vector>

std::vector<std::array<double, 3>>
makeStructuredCloud( double Lx, double Ly, double Lz, int nx, int ny, int nz,
                     double ox = 0., double oy = 0., double oz = 0. )
{
    std::vector<std::array<double, 3>> cloud( nx * ny * nz );
    std::function<int( int, int, int )> ind = [nx, ny]( int i, int j, int k ) {
        return i + j * nx + k * ( nx * ny );
    };
    double x, y, z;
    for ( int i = 0; i < nx; ++i )
        for ( int j = 0; j < ny; ++j )
            for ( int k = 0; k < nz; ++k )
            {
                x = ox + i * Lx / nx;
                y = oy + j * Ly / ny;
                z = oz + k * Lz / nz;
                cloud[ind( i, j, k )] = {{x, y, z}};
            }
    return cloud;
}

std::vector<std::array<double, 3>>
makeRandomCloud( double Lx, double Ly, double Lz, int n, double seed = 0. )
{
    std::vector<std::array<double, 3>> cloud( n );
    std::default_random_engine generator( seed );
    std::uniform_real_distribution<double> distributionx( 0.0, Lx );
    std::uniform_real_distribution<double> distributiony( 0.0, Ly );
    std::uniform_real_distribution<double> distributionz( 0.0, Lz );
    for ( int i = 0; i < n; ++i )
    {
        double x = distributionx( generator );
        double y = distributiony( generator );
        double z = distributionz( generator );
        cloud[i] = {{x, y, z}};
    }
    return cloud;
}

template <typename DeviceType>
void copyPointsFromCloud( std::vector<std::array<double, 3>> const &cloud,
                          Kokkos::View<double **, DeviceType> &points )
{
    int const n_points = cloud.size();
    int const spatial_dim = 3;
    Kokkos::realloc( points, n_points, spatial_dim );
    auto points_host = Kokkos::create_mirror_view( points );
    for ( int i = 0; i < n_points; ++i )
        for ( int d = 0; d < spatial_dim; ++d )
            points_host( i, d ) = cloud[i][d];
    Kokkos::deep_copy( points, points_host );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NearestNeighborOperator, unique_source_point,
                                   DeviceType )
{
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int const comm_size = comm->getSize();
    int const comm_rank = comm->getRank();

    int const space_dim = 3;

    // Build structured cloud of points for the source and random cloud for the
    // target.
    Kokkos::View<double **, DeviceType> source_points( "source" );

    Kokkos::View<double **, DeviceType> target_points( "target", 1, space_dim );

    TEST_THROW( DataTransferKit::NearestNeighborOperator<DeviceType>(
                    comm, source_points, target_points ),
                DataTransferKit::DataTransferKitException );

    if ( comm_rank == 0 )
    {
        Kokkos::resize( source_points, 1, space_dim );
        auto source_points_host = Kokkos::create_mirror_view( source_points );
        for ( int d = 0; d < space_dim; ++d )
            source_points_host( 0, d ) = (double)comm_size;
        Kokkos::deep_copy( source_points, source_points_host );
    }

    auto target_points_host = Kokkos::create_mirror_view( target_points );
    for ( int d = 0; d < space_dim; ++d )
        target_points_host( 0, d ) = (double)comm_rank;
    Kokkos::deep_copy( target_points, target_points_host );

    // Shameless hack to help the distributed tree with the nearest neighbor
    // search.
    auto const epsilon_default =
        DataTransferKit::DistributedSearchTreeImpl<DeviceType>::epsilon;
    DataTransferKit::DistributedSearchTreeImpl<DeviceType>::epsilon =
        (double)comm_size;

    DataTransferKit::NearestNeighborOperator<DeviceType> nnop(
        comm, source_points, target_points );

    Kokkos::View<double *, DeviceType> source_values( "in" );
    Kokkos::View<double *, DeviceType> target_values( "out" );

    // violate pre condition of apply
    TEST_THROW( nnop.apply( source_values, target_values ),
                DataTransferKit::DataTransferKitException );

    Kokkos::realloc( target_values, target_points.extent( 0 ) );
    Kokkos::realloc( source_values, source_points.extent( 0 ) );
    if ( comm_rank == 0 )
    {
        auto source_values_host = Kokkos::create_mirror_view( source_values );
        source_values_host( 0 ) = 255.;
        Kokkos::deep_copy( source_values, source_values_host );
    }

    nnop.apply( source_values, target_values );

    auto target_values_host = Kokkos::create_mirror_view( target_values );
    Kokkos::deep_copy( target_values_host, target_values );
    std::vector<double> target_values_ref = {255.};
    TEST_COMPARE_ARRAYS( target_values_host, target_values_ref );

    // Reset the static variable to its original value to avoid interfering
    // with other unit tests.
    DataTransferKit::DistributedSearchTreeImpl<DeviceType>::epsilon =
        epsilon_default;
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NearestNeighborOperator, structured_clouds,
                                   DeviceType )
{
    // The source is a structured cloud. The target is the same cloud but
    // distributed differently among the processors.
    Teuchos::RCP<Teuchos::Comm<int> const> comm =
        Teuchos::DefaultComm<int>::getComm();
    unsigned int const comm_size = comm->getSize();
    unsigned int const comm_rank = comm->getRank();

    // Build the structured cloud of points for the source and the target.
    double const Lx = 2.;
    double const Ly = 3.;
    double const Lz = 5.;
    unsigned int const nx = 7;
    unsigned int const ny = 11;
    unsigned int const nz = 13;
    double const source_offset_x = comm_rank * Lx;
    double const source_offset_y = comm_rank * Ly;
    double const source_offset_z = comm_rank * Lz;

    Kokkos::View<double **, DeviceType> source_points( "source_points" );
    copyPointsFromCloud<DeviceType>(
        makeStructuredCloud( Lx, Ly, Lz, nx, ny, nz, source_offset_x,
                             source_offset_y, source_offset_z ),
        source_points );

    double const target_offset_x = ( ( comm_rank + 1 ) % comm_size ) * Lx;
    double const target_offset_y = ( ( comm_rank + 1 ) % comm_size ) * Ly;
    double const target_offset_z = ( ( comm_rank + 1 ) % comm_size ) * Lz;

    Kokkos::View<double **, DeviceType> target_points( "target_points" );
    copyPointsFromCloud<DeviceType>(
        makeStructuredCloud( Lx, Ly, Lz, nx, ny, nz, target_offset_x,
                             target_offset_y, target_offset_z ),
        target_points );

    unsigned int const n_points = source_points.extent( 0 );
    Kokkos::View<double *, DeviceType> target_values( "target_values",
                                                      n_points );

    DataTransferKit::NearestNeighborOperator<DeviceType> nnop(
        comm, source_points, target_points );

    Kokkos::View<double *, DeviceType> source_values( "source_values",
                                                      n_points );
    Kokkos::deep_copy( source_values,
                       Kokkos::subview( source_points, Kokkos::ALL, 0 ) );

    nnop.apply( source_values, target_values );

    // Check results
    auto target_values_host = Kokkos::create_mirror_view( target_values );
    Kokkos::deep_copy( target_values_host, target_values );
    auto target_points_host = Kokkos::create_mirror_view( target_points );
    Kokkos::deep_copy( target_points_host, target_points );
    for ( unsigned int i = 0; i < n_points; ++i )
        TEST_FLOATING_EQUALITY( target_values_host( i ),
                                target_points_host( i, 0 ), 1e-14 );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( NearestNeighborOperator, mixed_clouds,
                                   DeviceType )
{
    // The source is a structured cloud. The target is a random cloud.
    Teuchos::RCP<Teuchos::Comm<int> const> comm =
        Teuchos::DefaultComm<int>::getComm();
    unsigned int const comm_size = comm->getSize();
    unsigned int const comm_rank = comm->getRank();

    // Build the structured cloud of points for the source.
    unsigned int constexpr spacedim = 3;
    double const Lx = 17.;
    double const Ly = 19.;
    double const Lz = 23.;
    unsigned int const nx = 29;
    unsigned int const ny = 31;
    unsigned int const nz = 37;
    double const source_offset_x = comm_rank * Lx;
    double const source_offset_y = 0;
    double const source_offset_z = 0;

    std::vector<std::array<double, spacedim>> structured_cloud =
        makeStructuredCloud( Lx, Ly, Lz, nx, ny, nz, source_offset_x,
                             source_offset_y, source_offset_z );
    Kokkos::View<double **, DeviceType> source_points( "source_points" );
    copyPointsFromCloud<DeviceType>( structured_cloud, source_points );

    // Build the random cloud of points for the target.
    unsigned int const n_target_points = 41;
    std::vector<std::array<double, spacedim>> random_cloud =
        makeRandomCloud( comm_size * Lx, Ly, Lz, n_target_points, comm_rank );

    Kokkos::View<double **, DeviceType> target_points( "target_points" );
    copyPointsFromCloud<DeviceType>( random_cloud, target_points );

    Kokkos::View<double *, DeviceType> target_values( "target_values",
                                                      n_target_points );

    // Approximate nearest neighbor search may fail in some situations:
    //
    //    X     X     X     X     X  source points
    //    ------------>     <------
    //    rank 0            rank 1
    //
    //                   ^
    //                   target point that falls in the gap does not overlap
    //                   with local trees so it won't be able to find its
    //                   neighbors
    auto const epsilon_default =
        DataTransferKit::DistributedSearchTreeImpl<DeviceType>::epsilon;
    TEST_THROW( DataTransferKit::NearestNeighborOperator<DeviceType>(
                    comm, source_points, target_points ),
                DataTransferKit::DataTransferKitException );

    // Determine appropriate tolerance for the approximate nearest neighbor
    // search on the distributed tree.
    DataTransferKit::DistributedSearchTreeImpl<DeviceType>::epsilon =
        std::max( {Lx / nx, Ly / ny, Lz / nz} );

    // This time we do not get the exception when we call the constructor.
    DataTransferKit::NearestNeighborOperator<DeviceType> nnop(
        comm, source_points, target_points );
    DataTransferKit::DistributedSearchTreeImpl<DeviceType>::epsilon =
        epsilon_default;

    unsigned int const n_points = source_points.extent( 0 );
    Kokkos::View<double *, DeviceType> source_values( "source_values",
                                                      n_points );
    Kokkos::deep_copy( source_values,
                       Kokkos::subview( source_points, Kokkos::ALL, 0 ) );

    nnop.apply( source_values, target_values );

    // Check results
    auto target_values_host = Kokkos::create_mirror_view( target_values );
    Kokkos::deep_copy( target_values_host, target_values );
    for ( unsigned int i = 0; i < n_target_points; ++i )
    {
        double ref_value = round( target_points( i, 0 ) / Lx * nx ) * Lx / nx;
        if ( ref_value == Lx * comm_size )
            ref_value -= Lx / nx;
        TEST_FLOATING_EQUALITY( target_values_host( i, 0 ), ref_value, 1e-14 );
    }
}

// Include the test macros.
#include "DataTransferKitMeshfree_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    using DeviceType##NODE = typename NODE::device_type;                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                      \
        NearestNeighborOperator, unique_source_point, DeviceType##NODE )       \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                      \
        NearestNeighborOperator, structured_clouds, DeviceType##NODE )         \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( NearestNeighborOperator,             \
                                          mixed_clouds, DeviceType##NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
