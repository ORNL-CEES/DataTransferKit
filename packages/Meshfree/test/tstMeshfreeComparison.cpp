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

#include <random>

constexpr int DIM = 3;

struct MeshFreeRunOptions
{
    MeshFreeRunOptions() = default;
    int knn = 0;
    double search_radius = 0.;
};

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
        grid_points.reserve( n_points[0] * n_points[1] * n_points[2] );

        std::array<double, DIM> h = {
            ( upper_corner[0] - lower_corner[0] ) / ( n_points[0] + 1 ),
            ( upper_corner[1] - lower_corner[1] ) / ( n_points[1] + 1 ),
            ( upper_corner[2] - lower_corner[2] ) / ( n_points[2] + 1 )};

        std::uniform_real_distribution<double> distribution( -0.5, 0.5 );
        std::default_random_engine generator;
        auto random = [&distribution, &generator]() {
            return 0.; // distribution( generator );
        };

        for ( int i = 0; i < n_points[0]; ++i )
            for ( int j = 0; j < n_points[1]; ++j )
                for ( int k = 0; k < n_points[2]; ++k )
                {
                    std::array<double, DIM> point{
                        lower_corner[0] + h[0] * ( i + 1 + random() ),
                        lower_corner[1] + h[1] * ( j + 1 + random() ),
                        lower_corner[2] + h[2] * ( k + 1 + random() )};
                    grid_points.push_back( point );
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
void testOperator( int source_points_per_dim, int target_points_per_dim,
                   MeshFreeRunOptions const &run_options )
{
    using DeviceType = typename Operator::device_type;

    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );

    int comm_size;
    MPI_Comm_size( comm, &comm_size );

    int l = std::ceil( std::cbrt( comm_size ) );

    // compute the offset necessary to fit local boxes in the unit cube.
    double offset_x = ( comm_rank % l ) * 1. / l;
    double offset_y = ( ( comm_rank / l ) % l ) * 1. / l;
    double offset_z = ( ( comm_rank / ( l * l ) ) % l ) * 1. / l;

    std::array<int, DIM> n_source_points_grid = {
        source_points_per_dim, source_points_per_dim, source_points_per_dim};
    std::array<double, DIM> lower_corner = {offset_x, offset_y, offset_z};
    std::array<double, DIM> upper_corner = {
        offset_x + 1. / l, offset_y + 1. / l, offset_z + 1. / l};
    auto source_points_arr = Helper<DeviceType>::makeGridPoints(
        n_source_points_grid, lower_corner, upper_corner );

    std::array<int, DIM> n_target_points_grid = {
        target_points_per_dim, target_points_per_dim, target_points_per_dim};
    auto target_points_arr = Helper<DeviceType>::makeGridPoints(
        n_target_points_grid, lower_corner, upper_corner );

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
               3 * p[1] * p[2] + 4 * p[2] * p[2] /*-
               std::sin( p[0] ) * std::sin( p[1] ) * std::sin( p[2] )*/
            ;
    };

    for ( unsigned int i = 0; i < n_source_points; ++i )
        source_values_arr[i] = f( source_points_arr[i] );
    for ( unsigned int i = 0; i < n_target_points; ++i )
        target_values_ref[i] = f( target_points_arr[i] );

    auto source_points = Helper<DeviceType>::makePoints( source_points_arr );
    auto source_values = Helper<DeviceType>::makeValues( source_values_arr );
    auto target_points = Helper<DeviceType>::makePoints( target_points_arr );
    auto target_values = Helper<DeviceType>::makeValues( target_values_arr );

    const unsigned int n_constructor_iterations = 10;
    std::unique_ptr<Operator> op_ptr;
    MPI_Barrier( comm );
    Kokkos::fence();
    auto start_setup = std::chrono::high_resolution_clock::now();
    for ( unsigned int i = 0; i < n_constructor_iterations; ++i )
    {
        op_ptr = std::make_unique<Operator>( comm, source_points, target_points,
                                             run_options.knn,
                                             run_options.search_radius );
        MPI_Barrier( comm );
        Kokkos::fence();
    }
    auto end_setup = std::chrono::high_resolution_clock::now();

    const unsigned int n_apply_iterations = 10;
    MPI_Barrier( comm );
    Kokkos::fence();
    auto start_apply = std::chrono::high_resolution_clock::now();
    for ( unsigned int i = 0; i < n_apply_iterations; ++i )
    {
        op_ptr->apply( source_values, target_values );
        MPI_Barrier( comm );
        Kokkos::fence();
    }
    auto end_apply = std::chrono::high_resolution_clock::now();

    auto target_values_host = Kokkos::create_mirror_view( target_values );
    Kokkos::deep_copy( target_values_host, target_values );

    double total_error = 0.;
    for ( unsigned int i = 0; i < target_values_host.extent( 0 ); ++i )
        total_error +=
            std::abs( target_values_host( i ) - target_values_ref[i] );
    total_error /= target_values_host.extent( 0 );
    double global_error = 0.;
    MPI_Allreduce( &total_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, comm );

    if ( comm_rank == 0 )
    {
        printf( "(%d,%d) knn: %d, radius: %e, setup: %lu ms, apply: %lu ms, "
                "error: %e\n",
                source_points_per_dim, target_points_per_dim, run_options.knn,
                run_options.search_radius,
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    end_setup - start_setup )
                        .count() /
                    n_constructor_iterations,
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    end_apply - start_apply )
                        .count() /
                    n_apply_iterations,
                global_error / comm_size );
    }
}

int main( int argc, char *argv[] )
{
    MPI_Init( &argc, &argv );
    Kokkos::initialize( argc, argv );

    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );

    using NODE = Kokkos::Serial;
    using Wendland = DataTransferKit::Wendland<0>;
    // using Wendland = DataTransferKit::Wendland<2>;
    // using Wendland = DataTransferKit::Wendland<6>;

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
        for ( unsigned int n_refinements = 0; n_refinements < 3;
              ++n_refinements )
        {
            MPI_Barrier( comm );
            if ( comm_rank == 0 )
                std::cout << "MLS 0    ";
            MeshFreeRunOptions run_options_mls_0{Constant3::size, 0.};
            testOperator<DataTransferKit::MovingLeastSquaresOperator<
                typename NODE::device_type, Wendland, Constant3>>(
                n_source_points, n_source_points / 2, run_options_mls_0 );

            MPI_Barrier( comm );

            if ( comm_rank == 0 )
                std::cout << "MLS 1    ";
            MeshFreeRunOptions run_options_mls_1{Linear3::size, 0.};
            testOperator<DataTransferKit::MovingLeastSquaresOperator<
                typename NODE::device_type, Wendland, Linear3>>(
                n_source_points, n_source_points / 2, run_options_mls_1 );

            MPI_Barrier( comm );

            if ( comm_rank == 0 )
                std::cout << "MLS 2    ";
            MeshFreeRunOptions run_options_mls_2{Quadratic3::size, 0.};
            testOperator<DataTransferKit::MovingLeastSquaresOperator<
                typename NODE::device_type, Wendland, Quadratic3>>(
                n_source_points, n_source_points / 2, run_options_mls_2 );

            MPI_Barrier( comm );

            if ( comm_rank == 0 )
                std::cout << "Spline 1 ";
            MeshFreeRunOptions run_options_spline_1{0, 1. / n_source_points};
            testOperator<DataTransferKit::SplineOperator<
                typename NODE::device_type, Wendland, Linear3>>(
                n_source_points, n_source_points / 2, run_options_spline_1 );

            n_source_points <<= 1;
        }
    }
    {
        int n_target_points = 2;
        for ( unsigned int n_refinements = 0; n_refinements < 3;
              ++n_refinements )
        {
            MPI_Barrier( comm );

            if ( comm_rank == 0 )
                std::cout << "MLS 0    ";
            MeshFreeRunOptions run_options_mls_0{Constant3::size, 0.};
            testOperator<DataTransferKit::MovingLeastSquaresOperator<
                typename NODE::device_type, Wendland, Constant3>>(
                n_target_points / 2, n_target_points, run_options_mls_0 );

            MPI_Barrier( comm );

            if ( comm_rank == 0 )
                std::cout << "MLS 1    ";
            MeshFreeRunOptions run_options_mls_1{Linear3::size, 0.};
            testOperator<DataTransferKit::MovingLeastSquaresOperator<
                typename NODE::device_type, Wendland, Linear3>>(
                n_target_points / 2, n_target_points, run_options_mls_1 );

            MPI_Barrier( comm );

            if ( comm_rank == 0 )
                std::cout << "MLS 2    ";
            MeshFreeRunOptions run_options_mls_2{Quadratic3::size, 0.};
            testOperator<DataTransferKit::MovingLeastSquaresOperator<
                typename NODE::device_type, Wendland, Quadratic3>>(
                n_target_points / 2, n_target_points, run_options_mls_2 );

            MPI_Barrier( comm );

            if ( comm_rank == 0 )
                std::cout << "Spline 1 ";
            MeshFreeRunOptions run_options_spline_1{0, 4. / n_target_points};
            testOperator<DataTransferKit::SplineOperator<
                typename NODE::device_type, Wendland, Linear3>>(
                n_target_points / 2, n_target_points, run_options_spline_1 );

            n_target_points <<= 1;
        }
    }

    if ( comm_rank == 0 )
        std::cout << "End Result: TEST PASSED" << std::endl;

    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}
