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

#include <Teuchos_UnitTestHarness.hpp>

#include <DTK_DBC.hpp> // DataTransferKitException
#include <DTK_MovingLeastSquaresOperator_decl.hpp>
#include <DTK_MovingLeastSquaresOperator_def.hpp>
#include <DTK_SplineOperator_decl.hpp>
#include <DTK_SplineOperator_def.hpp>
#include <Kokkos_Core.hpp>

#include <array>
#include <cmath>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

int constexpr DIM = 3;

struct MLS
{
};
struct Spline
{
};

template <typename DeviceType>
struct Helper
{

    static void
    makeSourceTargetPoints( std::vector<std::array<double, DIM>> &source_points,
                            std::vector<std::array<double, DIM>> &target_points,
                            int n_source_points_in_radius, double radius,
                            int rank )
    {
        const int n_target_points = target_points.size();

        std::default_random_engine g;
        std::uniform_real_distribution<double> rx( -radius, radius );

        double const offset = n_target_points * radius * 10. * rank;
        for ( int i = 0; i < n_target_points; i++ )
        {
            // c is the center of our small universe
            std::array<double, DIM> c = {{i * ceil( 4 * radius ) + offset,
                                          i * ceil( 4 * radius ) + offset,
                                          i * ceil( 4 * radius ) + offset}};

            // Construct a target point that is within [-0.5*radius, 0.5*radius]
            // of the center in each dimension
            target_points[i] = {{c[0] + 0.5 * rx( g ), c[1] + 0.5 * rx( g ),
                                 c[2] + 0.5 * rx( g )}};

            // Construct source points that are within [-radius, radius] of
            // the center in each dimension. This would mean that the distance
            // between any of these source points and target point is no more
            // than 2*radius.
            for ( int k = 0; k < n_source_points_in_radius; k++ )
            {
                auto ind = i * n_source_points_in_radius + k;

                source_points[ind] = {
                    {c[0] + rx( g ), c[1] + rx( g ), c[2] + rx( g )}};
            }
        }
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

    static std::vector<std::array<double, DIM>>
    makeGridPoints( std::array<int, DIM> const &n_points,
                    std::array<double, DIM> const &offset )
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
                        offset[0] + i, offset[1] + j, offset[2] + k};
                    grid_points.push_back( point );
                }
            }
        }

        return grid_points;
    }
};

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MeshfreeOperator, same_npoints_and_basis,
                                   OperatorType, Operator )
{
    // Test the situation when the number of found source points is the same as
    // the basis size. This greatly simplifies the situation as the equation
    //   P^T PHI P a = P^T PHI f
    // becomes
    //   P a = f
    // meaning that as long as the points are "independent", the manifold
    // corresponding to the order of the basis is going to be approximated
    // exactly. In other words, the result of the approximation is going to be
    // the exact value of the function. For instance, for Linear basis
    // (1,x,y,z) in 3D, any linear function would be approximated exactly
    // independent of the chosen radial basis function.
    //
    // We perform the test in a "batched" version meaning we simultaneously
    // test multiple target points. The test is constructed in such a way that
    // the number of source points found within the specified distance of each
    // target point is exactly the size of the desired basis.
    //
    // NOTE: right now, the DIM = 3 is hardcoded.
    using namespace DataTransferKit;

    using DeviceType = typename Operator::device_type;
    using PolynomialBasis = typename Operator::polynomial_basis;

    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );

    const int n_target_points = 10;
    const double radius = 1.0;

    const int n_source_points_in_radius = PolynomialBasis::size;

    const int n_source_points = n_target_points * n_source_points_in_radius;

    std::vector<std::array<double, DIM>> source_points_arr( n_source_points );
    std::vector<std::array<double, DIM>> target_points_arr( n_target_points );

    std::vector<double> source_values_arr( n_source_points );
    std::vector<double> target_values_arr( n_target_points );
    std::vector<double> target_values_ref( n_target_points );

    Helper<DeviceType>::makeSourceTargetPoints(
        source_points_arr, target_points_arr, n_source_points_in_radius,
        0.5 * radius, comm_rank );

    // Arbitrary function of the specified order
    std::function<double( std::array<double, DIM> )> f;
    switch ( PolynomialBasis::size )
    {
    case 1: // constant
        f = []( std::array<double, DIM> ) -> double { return 3.0; };
        break;
    case 4: // linear
        f = []( std::array<double, DIM> p ) -> double {
            return 4 + 2 * p[0] + 3 * p[1] - 2 * p[2];
        };
        break;
    case 10: // quadratic
        f = []( std::array<double, DIM> p ) -> double {
            return 2 + 3 * p[0] - 5 * p[1] + 2 * p[2] + 3 * p[0] * p[0] +
                   4 * p[0] * p[1] - 2 * p[0] * p[2] + p[1] * p[1] -
                   3 * p[1] * p[2] + 4 * p[2] * p[2];
        };
        break;
    default:
        throw;
    };

    for ( int i = 0; i < n_source_points; i++ )
        source_values_arr[i] = f( source_points_arr[i] );
    for ( int i = 0; i < n_target_points; i++ )
        target_values_ref[i] = f( target_points_arr[i] );

    auto source_points = Helper<DeviceType>::makePoints( source_points_arr );
    auto source_values = Helper<DeviceType>::makeValues( source_values_arr );
    auto target_points = Helper<DeviceType>::makePoints( target_points_arr );
    auto target_values = Helper<DeviceType>::makeValues( target_values_arr );

    Operator op( comm, source_points, target_points );

    op.apply( source_values, target_values );

    double eps = 0.0;
    if ( std::is_same<OperatorType, MLS>{} )
        eps = 1e-7;
    else if ( std::is_same<OperatorType, Spline>{} )
        eps = 2e-7;

    auto target_values_host = Kokkos::create_mirror_view( target_values );
    Kokkos::deep_copy( target_values_host, target_values );
    TEST_COMPARE_FLOATING_ARRAYS( target_values_host, target_values_ref, eps );
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MeshfreeOperator, grid, OperatorType,
                                   Operator )
{
    using namespace DataTransferKit;

    using DeviceType = typename Operator::device_type;
    using PolynomialBasis = typename Operator::polynomial_basis;

    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );

    std::array<int, DIM> n_source_points_grid = {40, 40, 1};
    std::array<double, DIM> offset = {0., 0., static_cast<double>( comm_rank )};
    auto source_points_arr =
        Helper<DeviceType>::makeGridPoints( n_source_points_grid, offset );

    std::array<int, DIM> n_target_points_grid = {1, 1, 1};
    offset = {19, 19., static_cast<double>( comm_rank )};
    auto target_points_arr =
        Helper<DeviceType>::makeGridPoints( n_target_points_grid, offset );

    unsigned int const n_source_points = source_points_arr.size();
    unsigned int const n_target_points = target_points_arr.size();
    std::vector<double> source_values_arr( n_source_points );
    std::vector<double> target_values_arr( n_target_points );
    std::vector<double> target_values_ref( n_target_points );

    // Arbitrary function of the specified order
    std::function<double( std::array<double, DIM> )> f;
    switch ( PolynomialBasis::size )
    {
    case 1: // constant
        f = []( std::array<double, DIM> ) -> double { return 3.0; };
        break;
    case 4: // linear
        f = []( std::array<double, DIM> p ) -> double {
            return 4 + 2 * p[0] + 3 * p[1] - 2 * p[2];
        };
        break;
    case 10: // quadratic
        f = []( std::array<double, DIM> p ) -> double {
            return 2 + 3 * p[0] - 5 * p[1] + 2 * p[2] + 3 * p[0] * p[0] +
                   4 * p[0] * p[1] - 2 * p[0] * p[2] + p[1] * p[1] -
                   3 * p[1] * p[2] + 4 * p[2] * p[2];
        };
        break;
    default:
        throw;
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

    double eps = 0.0;
    if ( std::is_same<OperatorType, MLS>{} )
        eps = 1e-14;
    else if ( std::is_same<OperatorType, Spline>{} )
        eps = 1e-9;

    auto target_values_host = Kokkos::create_mirror_view( target_values );
    Kokkos::deep_copy( target_values_host, target_values );
    TEST_COMPARE_FLOATING_ARRAYS( target_values_host, target_values_ref, eps );
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MeshfreeOperator, line, OperatorType,
                                   Operator )
{
    using namespace DataTransferKit;

    using DeviceType = typename Operator::device_type;
    using PolynomialBasis = typename Operator::polynomial_basis;

    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );

    std::array<int, DIM> n_source_points_grid = {40, 1, 1};
    std::array<double, DIM> offset = {0., 0., static_cast<double>( comm_rank )};
    auto source_points_arr =
        Helper<DeviceType>::makeGridPoints( n_source_points_grid, offset );

    std::array<int, DIM> n_target_points_grid = {39, 1, 1};
    offset = {0.5, 0., static_cast<double>( comm_rank )};
    auto target_points_arr =
        Helper<DeviceType>::makeGridPoints( n_target_points_grid, offset );

    unsigned int const n_source_points = source_points_arr.size();
    unsigned int const n_target_points = target_points_arr.size();
    std::vector<double> source_values_arr( n_source_points );
    std::vector<double> target_values_arr( n_target_points );
    std::vector<double> target_values_ref( n_target_points );

    // Arbitrary function of the specified order
    std::function<double( std::array<double, DIM> )> f;
    switch ( PolynomialBasis::size )
    {
    case 1: // constant
        f = []( std::array<double, DIM> ) -> double { return 3.0; };
        break;
    case 4: // linear
        f = []( std::array<double, DIM> p ) -> double {
            return 4 + 2 * p[0] + 3 * p[1] - 2 * p[2];
        };
        break;
    case 10: // quadratic
        f = []( std::array<double, DIM> p ) -> double {
            return 2 + 3 * p[0] - 5 * p[1] + 2 * p[2] + 3 * p[0] * p[0] +
                   4 * p[0] * p[1] - 2 * p[0] * p[2] + p[1] * p[1] -
                   3 * p[1] * p[2] + 4 * p[2] * p[2];
        };
        break;
    default:
        throw;
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

    double eps = 0.0;
    if ( std::is_same<OperatorType, MLS>{} )
        eps = 1e-14;
    else if ( std::is_same<OperatorType, Spline>{} )
        eps = 2e-9;

    auto target_values_host = Kokkos::create_mirror_view( target_values );
    Kokkos::deep_copy( target_values_host, target_values );
    TEST_COMPARE_FLOATING_ARRAYS( target_values_host, target_values_ref, eps );
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MeshfreeOperator, single_point_in_radius,
                                   OperatorType, Operator )
{
    using namespace DataTransferKit;

    using DeviceType = typename Operator::device_type;
    using PolynomialBasis = typename Operator::polynomial_basis;

    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_rank;
    MPI_Comm_rank( comm, &comm_rank );

    const int n_target_points = 10;
    const double radius = 1.0;

    const int n_source_points_in_radius = 1;

    const int n_source_points = n_target_points * n_source_points_in_radius;

    std::vector<std::array<double, DIM>> source_points_arr( n_source_points );
    std::vector<std::array<double, DIM>> target_points_arr( n_target_points );

    std::vector<double> source_values_arr( n_source_points );
    std::vector<double> target_values_arr( n_target_points );

    Helper<DeviceType>::makeSourceTargetPoints(
        source_points_arr, target_points_arr, n_source_points_in_radius,
        0.5 * radius, comm_rank );

    // Arbitrary function of the specified order
    std::function<double( std::array<double, DIM> )> f;
    switch ( PolynomialBasis::size )
    {
    case 1: // constant
        f = []( std::array<double, DIM> ) -> double { return 3.0; };
        break;
    case 4: // linear
        f = []( std::array<double, DIM> p ) -> double {
            return 4 + 2 * p[0] + 3 * p[1] - 2 * p[2];
        };
        break;
    case 10: // quadratic
        f = []( std::array<double, DIM> p ) -> double {
            return 2 + 3 * p[0] - 5 * p[1] + 2 * p[2] + 3 * p[0] * p[0] +
                   4 * p[0] * p[1] - 2 * p[0] * p[2] + p[1] * p[1] -
                   3 * p[1] * p[2] + 4 * p[2] * p[2];
        };
        break;
    default:
        throw;
    };

    for ( int i = 0; i < n_source_points; i++ )
        source_values_arr[i] = f( source_points_arr[i] );

    auto source_points = Helper<DeviceType>::makePoints( source_points_arr );
    auto source_values = Helper<DeviceType>::makeValues( source_values_arr );
    auto target_points = Helper<DeviceType>::makePoints( target_points_arr );
    auto target_values = Helper<DeviceType>::makeValues( target_values_arr );

    Operator op( comm, source_points, target_points );

    op.apply( source_values, target_values );

    // The approximation is really poor so we just check that the solution is
    // finite.
    auto target_values_host = Kokkos::create_mirror_view( target_values );
    Kokkos::deep_copy( target_values_host, target_values );
    for ( unsigned int i = 0; i < target_values_host.extent( 0 ); ++i )
        TEST_ASSERT( std::isfinite( target_values_host[i] ) );
}

// Include the test macros.
#include "DataTransferKit_ETIHelperMacros.h"

using Wendland0 = DataTransferKit::Wendland<0>;
using Wendland2 = DataTransferKit::Wendland<2>;
using Wendland6 = DataTransferKit::Wendland<6>;
using Constant3 =
    DataTransferKit::MultivariatePolynomialBasis<DataTransferKit::Constant, 3>;
using Linear3 =
    DataTransferKit::MultivariatePolynomialBasis<DataTransferKit::Linear, 3>;
using Quadratic3 =
    DataTransferKit::MultivariatePolynomialBasis<DataTransferKit::Quadratic, 3>;

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    using MLS_Wendland0_Constant3_##NODE =                                     \
        DataTransferKit::MovingLeastSquaresOperator<                           \
            typename NODE::device_type, Wendland0, Constant3>;                 \
    using MLS_Wendland0_Linear3_##NODE =                                       \
        DataTransferKit::MovingLeastSquaresOperator<                           \
            typename NODE::device_type, Wendland0, Linear3>;                   \
    using MLS_Wendland0_Quadratic3_##NODE =                                    \
        DataTransferKit::MovingLeastSquaresOperator<                           \
            typename NODE::device_type, Wendland0, Quadratic3>;                \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator,                    \
                                          same_npoints_and_basis, MLS,         \
                                          MLS_Wendland0_Constant3_##NODE )     \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator,                    \
                                          same_npoints_and_basis, MLS,         \
                                          MLS_Wendland0_Linear3_##NODE )       \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator,                    \
                                          same_npoints_and_basis, MLS,         \
                                          MLS_Wendland0_Quadratic3_##NODE )    \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator, line, MLS,         \
                                          MLS_Wendland0_Constant3_##NODE )     \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator, line, MLS,         \
                                          MLS_Wendland0_Linear3_##NODE )       \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator, line, MLS,         \
                                          MLS_Wendland0_Quadratic3_##NODE )    \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator, grid, MLS,         \
                                          MLS_Wendland0_Constant3_##NODE )     \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator, grid, MLS,         \
                                          MLS_Wendland0_Linear3_##NODE )       \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator, grid, MLS,         \
                                          MLS_Wendland0_Quadratic3_##NODE )    \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator,                    \
                                          single_point_in_radius, MLS,         \
                                          MLS_Wendland0_Constant3_##NODE )     \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator,                    \
                                          single_point_in_radius, MLS,         \
                                          MLS_Wendland0_Linear3_##NODE )       \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator,                    \
                                          single_point_in_radius, MLS,         \
                                          MLS_Wendland0_Quadratic3_##NODE )    \
    using Spline_Wendland0_Linear3_##NODE =                                    \
        DataTransferKit::SplineOperator<typename NODE::device_type, Wendland0, \
                                        Linear3>;                              \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator,                    \
                                          same_npoints_and_basis, Spline,      \
                                          Spline_Wendland0_Linear3_##NODE )    \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator, line, Spline,      \
                                          Spline_Wendland0_Linear3_##NODE )    \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator, grid, Spline,      \
                                          Spline_Wendland0_Linear3_##NODE )    \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MeshfreeOperator,                    \
                                          single_point_in_radius, Spline,      \
                                          Spline_Wendland0_Linear3_##NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
