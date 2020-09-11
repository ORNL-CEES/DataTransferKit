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

#include <DTK_MovingLeastSquaresOperator.hpp>
#include <DTK_SplineOperator.hpp>

constexpr int DIM = 3;

template <typename DeviceType>
struct Helper
{
    static std::vector<std::array<double, DIM>>
    makeGridPoints( std::array<int, DIM> const &n_points,
                    std::array<double, DIM> const &lower_corner,
                    std::array<double, DIM> const &upper_corner )
    {
        static_assert( DIM == 3, "Assume three dimensional geometry" );

        std::vector<std::array<double, DIM>> grid_points;

        for ( int i = 0; i < n_points[0]; ++i )
        {
            for ( int j = 0; j < n_points[1]; ++j )
            {
                for ( int k = 0; k < n_points[2]; ++k )
                {
                    std::array<double, DIM> point = {
                        lower_corner[0] +
                            ( upper_corner[0] - lower_corner[0] ) /
                                ( n_points[0] + 1 ) * ( i + 1 ),
                        lower_corner[1] +
                            ( upper_corner[1] - lower_corner[1] ) /
                                ( n_points[1] + 1 ) * ( j + 1 ),
                        lower_corner[2] +
                            ( upper_corner[2] - lower_corner[2] ) /
                                ( n_points[2] + 1 ) * ( k + 1 )};
                    grid_points.push_back( point );
                }
            }
        }

        return grid_points;
    }

    static Kokkos::View<DataTransferKit::Coordinate **, DeviceType>
    makePoints( std::vector<std::array<double, DIM>> const &in )
    {
        int const n = in.size();
        Kokkos::View<DataTransferKit::Coordinate **, DeviceType> out( "points",
                                                                      n, DIM );
        auto out_host = Kokkos::create_mirror_view( out );
        for ( int i = 0; i < n; ++i )
            for ( int j = 0; j < DIM; ++j )
                out_host( i, j ) = in[i][j];
        Kokkos::deep_copy( out, out_host );
        return out;
    }

    static Kokkos::View<double *, DeviceType>
    makeValues( std::vector<double> const &in )
    {
        int const n = in.size();
        Kokkos::View<double *, DeviceType> out( "points", n );
        auto out_host = Kokkos::create_mirror_view( out );
        for ( int i = 0; i < n; ++i )
            out_host( i ) = in[i];
        Kokkos::deep_copy( out, out_host );
        return out;
    }
};

template <typename Operator>
void testOperator( int source_points_per_dim, int target_points_per_dim )
{
    using DeviceType = typename Operator::device_type;

    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );

    std::array<int, DIM> n_source_points_grid = {
        source_points_per_dim, source_points_per_dim, source_points_per_dim};
    // FIXME MPI
    std::array<double, DIM> lower = {0., 0., static_cast<double>( comm_rank )};
    std::array<double, DIM> upper = {1., 1., 1.};
    auto source_points_arr = Helper<DeviceType>::makeGridPoints(
        n_source_points_grid, lower, upper );

    std::array<int, DIM> n_target_points_grid = {
        target_points_per_dim, target_points_per_dim, target_points_per_dim};
    auto target_points_arr = Helper<DeviceType>::makeGridPoints(
        n_target_points_grid, lower, upper );

    unsigned int const n_source_points = source_points_arr.size();
    unsigned int const n_target_points = target_points_arr.size();
    std::vector<double> source_values_arr( n_source_points );
    std::vector<double> target_values_arr( n_target_points );
    std::vector<double> target_values_ref( n_target_points );

    // Arbitrary function of the specified order
    std::function<double( std::array<double, DIM> )> f =
        []( std::array<double, DIM> p ) -> double {
        return 2 + 3 * p[0] - 5 * p[1] + 2 * p[2] + 3 * p[0] * p[0] +
               4 * p[0] * p[1] - 2 * p[0] * p[2] + p[1] * p[1] -
               3 * p[1] * p[2] + 4 * p[2] * p[2] -
               std::sin( p[0] ) * std::sin( p[1] ) * std::sin( p[2] );
    };

    for ( unsigned int i = 0; i < n_source_points; ++i )
        source_values_arr[i] = f( source_points_arr[i] );
    for ( unsigned int i = 0; i < n_target_points; ++i )
        target_values_ref[i] = f( target_points_arr[i] );

    auto source_points = Helper<DeviceType>::makePoints( source_points_arr );
    auto source_values = Helper<DeviceType>::makeValues( source_values_arr );
    auto target_points = Helper<DeviceType>::makePoints( target_points_arr );
    auto target_values = Helper<DeviceType>::makeValues( target_values_arr );

    Operator op( comm, source_points, target_points );

    op.apply( source_values, target_values );

    auto target_values_host = Kokkos::create_mirror_view( target_values );
    Kokkos::deep_copy( target_values_host, target_values );

    double total_error = 0.;
    for ( unsigned int i = 0; i < target_values_host.extent( 0 ); ++i )
        total_error =
            +std::abs( target_values_host( i ) - target_values_ref[i] );
    total_error /= target_values_host.extent( 0 );
    std::cout << "(" << source_points_per_dim << "," << target_points_per_dim
              << ") error: " << total_error << " ";
}

int main( int argc, char *argv[] )
{
    MPI_Init( &argc, &argv );
    Kokkos::initialize( argc, argv );

    using NODE = Kokkos::Serial;
    using Wendland0 = DataTransferKit::Wendland<0>;
    // using Wendland2 = DataTransferKit::Wendland<2>;
    // using Wendland6 = DataTransferKit::Wendland<6>;

    using Constant3 =
        DataTransferKit::MultivariatePolynomialBasis<DataTransferKit::Constant,
                                                     3>;
    using Linear3 =
        DataTransferKit::MultivariatePolynomialBasis<DataTransferKit::Linear,
                                                     3>;
    using Quadratic3 =
        DataTransferKit::MultivariatePolynomialBasis<DataTransferKit::Quadratic,
                                                     3>;

    {
        int n_source_points = 2;
        for ( unsigned int n_refinements = 0; n_refinements < 7;
              ++n_refinements )
        {
            std::cout << "MLS 0    ";
            auto t0 = std::chrono::high_resolution_clock::now();
            testOperator<DataTransferKit::MovingLeastSquaresOperator<
                typename NODE::device_type, Wendland0, Constant3>>(
                n_source_points, n_source_points / 2 );
            auto t1 = std::chrono::high_resolution_clock::now();
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
                             t1 - t0 )
                             .count()
                      << " ms" << std::endl;

            std::cout << "MLS 1    ";
            auto t2 = std::chrono::high_resolution_clock::now();
            testOperator<DataTransferKit::MovingLeastSquaresOperator<
                typename NODE::device_type, Wendland0, Linear3>>(
                n_source_points, n_source_points / 2 );
            auto t3 = std::chrono::high_resolution_clock::now();
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
                             t3 - t2 )
                             .count()
                      << " ms" << std::endl;

            std::cout << "MLS 2    ";
            auto t4 = std::chrono::high_resolution_clock::now();
            testOperator<DataTransferKit::MovingLeastSquaresOperator<
                typename NODE::device_type, Wendland0, Quadratic3>>(
                n_source_points, n_source_points / 2 );
            auto t5 = std::chrono::high_resolution_clock::now();
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
                             t5 - t4 )
                             .count()
                      << " ms" << std::endl;

            std::cout << "Spline 1 ";
            auto t6 = std::chrono::high_resolution_clock::now();
            testOperator<DataTransferKit::SplineOperator<
                typename NODE::device_type, Wendland0, Linear3>>(
                n_source_points, n_source_points / 2 );
            auto t7 = std::chrono::high_resolution_clock::now();
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
                             t7 - t6 )
                             .count()
                      << " ms" << std::endl;

            n_source_points <<= 1;
        }
    }
    {
        int n_target_points = 2;
        for ( unsigned int n_refinements = 0; n_refinements < 6;
              ++n_refinements )
        {
            std::cout << "MLS 0    ";
            auto t0 = std::chrono::high_resolution_clock::now();
            testOperator<DataTransferKit::MovingLeastSquaresOperator<
                typename NODE::device_type, Wendland0, Constant3>>(
                n_target_points / 2, n_target_points );
            auto t1 = std::chrono::high_resolution_clock::now();
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
                             t1 - t0 )
                             .count()
                      << " ms" << std::endl;

            std::cout << "MLS 1    ";
            auto t2 = std::chrono::high_resolution_clock::now();
            testOperator<DataTransferKit::MovingLeastSquaresOperator<
                typename NODE::device_type, Wendland0, Linear3>>(
                n_target_points / 2, n_target_points );
            auto t3 = std::chrono::high_resolution_clock::now();
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
                             t3 - t2 )
                             .count()
                      << " ms" << std::endl;

            std::cout << "MLS 2    ";
            auto t4 = std::chrono::high_resolution_clock::now();
            testOperator<DataTransferKit::MovingLeastSquaresOperator<
                typename NODE::device_type, Wendland0, Quadratic3>>(
                n_target_points / 2, n_target_points );
            auto t5 = std::chrono::high_resolution_clock::now();
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
                             t5 - t4 )
                             .count()
                      << " ms" << std::endl;

            std::cout << "Spline 1 ";
            auto t6 = std::chrono::high_resolution_clock::now();
            testOperator<DataTransferKit::SplineOperator<
                typename NODE::device_type, Wendland0, Linear3>>(
                n_target_points / 2, n_target_points );
            auto t7 = std::chrono::high_resolution_clock::now();
            std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
                             t7 - t6 )
                             .count()
                      << " ms" << std::endl;

            n_target_points <<= 1;
        }
    }

    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}
