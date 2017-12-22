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

#include <DTK_DBC.hpp> // DataTransferKitException
#include <DTK_MovingLeastSquaresOperator_decl.hpp>
#include <DTK_MovingLeastSquaresOperator_def.hpp>
#include <Kokkos_Core.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>

#include <array>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

template <typename DeviceType>
struct Helper
{
    static int const DIM = 3;

    static void
    makeSourceTargetPoints( std::vector<std::array<double, DIM>> &source_points,
                            std::vector<std::array<double, DIM>> &target_points,
                            int n_source_points_in_radius, double radius )
    {
        const int n_target_points = target_points.size();

        std::default_random_engine g;
        std::uniform_real_distribution<double> rx( -radius, radius );

        for ( int i = 0; i < n_target_points; i++ )
        {
            // c is the center of our small universe
            std::array<double, DIM> c = {{i * ceil( 4 * radius ),
                                          i * ceil( 4 * radius ),
                                          i * ceil( 4 * radius )}};

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

    static Kokkos::View<double **, DeviceType>
    makePoints( std::vector<std::array<double, DIM>> const &in )
    {
        int const n = in.size();
        Kokkos::View<double **, DeviceType> out( "points", n, DIM );
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

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MovingLeastSquaresOperator,
                                   same_npoints_and_basis, DeviceType,
                                   RadialBasisFunction, PolynomialBasis )
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

    auto comm = Teuchos::DefaultComm<int>::getComm();

    const auto DIM = Helper<DeviceType>::DIM;

    const int n_target_points = 10;
    const double radius = 1.0;
    const double R = sqrt( DIM ) * radius;

    const int size_polynomial_basis = PolynomialBasis::size();
    const int n_source_points_in_radius = size_polynomial_basis;

    const int n_source_points = n_target_points * n_source_points_in_radius;

    std::vector<std::array<double, DIM>> source_points_arr( n_source_points );
    std::vector<std::array<double, DIM>> target_points_arr( n_target_points );

    std::vector<double> source_values_arr( n_source_points );
    std::vector<double> target_values_arr( n_target_points );
    std::vector<double> target_values_ref( n_target_points );

    Helper<DeviceType>::makeSourceTargetPoints(
        source_points_arr, target_points_arr, n_source_points_in_radius,
        0.5 * radius );

    // Arbitrary function of the specified order
    std::function<double( std::array<double, DIM> )> f;
    switch ( size_polynomial_basis )
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

    Teuchos::ParameterList plist;
    plist.set( "radius", R );

    DataTransferKit::MovingLeastSquaresOperator<DeviceType, RadialBasisFunction,
                                                PolynomialBasis>
        mlsop( comm, source_points, target_points, plist );

    mlsop.apply( source_values, target_values );

    TEST_COMPARE_FLOATING_ARRAYS( target_values, target_values_ref, 1e-11 );
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MovingLeastSquaresOperator, constant_basis,
                                   DeviceType, RadialBasisFunction )
{
    // Test the situation when the number of found source points is greater
    // than the basis size. Using constant basis greatly simplifies the
    // situation as the equation
    //   P^T PHI P a = P^T PHI f
    // becomes (due to P = [1, ..., 1]^T
    //   (\sum_k phi_k) a = \sum_k phi_k f_k
    // meaning that the polynomial value a (which is equal to the target value
    // in this case) is the weighted average of source values.
    //
    // We perform the test in a "batched" version meaning we simultaneously
    // test multiple target points. The test is constructed in such a way that
    // the number of source points found within the specified distance of each
    // target point is exactly the specified number which is larger than 1.
    //
    // NOTE: right now, the DIM = 3 is hardcoded.
    using namespace DataTransferKit;

    using PolynomialBasis = MultivariatePolynomialBasis<Constant, 3>;

    auto comm = Teuchos::DefaultComm<int>::getComm();

    const auto DIM = Helper<DeviceType>::DIM;

    const int n_target_points = 13;
    const double radius = 1.0;
    const double R = sqrt( DIM ) * radius;

    const int n_source_points_in_radius = 7;

    const int n_source_points = n_target_points * n_source_points_in_radius;

    std::vector<std::array<double, DIM>> source_points_arr( n_source_points );
    std::vector<std::array<double, DIM>> target_points_arr( n_target_points );

    std::vector<double> source_values_arr( n_source_points );
    std::vector<double> target_values_arr( n_target_points );
    std::vector<double> target_values_ref( n_target_points );

    Helper<DeviceType>::makeSourceTargetPoints(
        source_points_arr, target_points_arr, n_source_points_in_radius,
        0.5 * radius );

    // Arbitrary function
    std::function<double( std::array<double, DIM> )> f =
        []( std::array<double, DIM> p ) -> double {
        return sin( p[0] ) - 2 * cos( p[1] ) + p[2];
    };

    for ( int i = 0; i < n_source_points; i++ )
        source_values_arr[i] = f( source_points_arr[i] );

    const auto rbf = RadialBasisFunction( R );
    for ( int i = 0; i < n_target_points; i++ )
    {
        double phi_sum = 0.0;
        double f_sum = 0.0;
        for ( int k = 0; k < n_source_points_in_radius; k++ )
        {
            auto ind = i * n_source_points_in_radius + k;

            double phi = rbf( Details::distance(
                Point{{source_points_arr[ind][0], source_points_arr[ind][1],
                       source_points_arr[ind][2]}},
                Point{{target_points_arr[i][0], target_points_arr[i][1],
                       target_points_arr[i][2]}} ) );
            phi_sum += phi;
            f_sum += phi * source_values_arr[ind];
        }

        target_values_ref[i] = f_sum / phi_sum;
    }

    auto source_points = Helper<DeviceType>::makePoints( source_points_arr );
    auto source_values = Helper<DeviceType>::makeValues( source_values_arr );
    auto target_points = Helper<DeviceType>::makePoints( target_points_arr );
    auto target_values = Helper<DeviceType>::makeValues( target_values_arr );

    Teuchos::ParameterList plist;
    plist.set( "radius", R );

    DataTransferKit::MovingLeastSquaresOperator<DeviceType, RadialBasisFunction,
                                                PolynomialBasis>
        mlsop( comm, source_points, target_points, plist );

    mlsop.apply( source_values, target_values );

    TEST_COMPARE_FLOATING_ARRAYS( target_values, target_values_ref, 1e-13 );
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MovingLeastSquaresOperator, underdetermined,
                                   DeviceType, RadialBasisFunction,
                                   PolynomialBasis )
{
    // Test the situation when the number of found source points is less than
    // the basis size. This is the toughest case as it makes the system
    //   P^T PHI P a = P^T PHI f
    // be of not full rank.
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

    auto comm = Teuchos::DefaultComm<int>::getComm();

    const auto DIM = Helper<DeviceType>::DIM;

    const int n_target_points = 10;
    const double radius = 1.0;
    const double R = sqrt( DIM ) * radius;

    const int size_polynomial_basis = PolynomialBasis::size();
    const int n_source_points_in_radius = 1;

    const int n_source_points = n_target_points * n_source_points_in_radius;

    std::vector<std::array<double, DIM>> source_points_arr( n_source_points );
    std::vector<std::array<double, DIM>> target_points_arr( n_target_points );

    std::vector<double> source_values_arr( n_source_points );
    std::vector<double> target_values_arr( n_target_points );
    std::vector<double> target_values_ref( n_target_points );

    Helper<DeviceType>::makeSourceTargetPoints(
        source_points_arr, target_points_arr, n_source_points_in_radius,
        0.5 * radius );

    // Arbitrary function of the specified order
    std::function<double( std::array<double, DIM> )> f;
    switch ( size_polynomial_basis )
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

    if ( n_source_points_in_radius == 1 )
    {
        for ( int i = 0; i < n_target_points; i++ )
            target_values_ref[i] = source_values_arr[i];
    }

    auto source_points = Helper<DeviceType>::makePoints( source_points_arr );
    auto source_values = Helper<DeviceType>::makeValues( source_values_arr );
    auto target_points = Helper<DeviceType>::makePoints( target_points_arr );
    auto target_values = Helper<DeviceType>::makeValues( target_values_arr );

    Teuchos::ParameterList plist;
    plist.set( "radius", R );

    // FIXME: we currently have no robust way to deal with underdetermined
    // systems. So, we just do the stupid thing and throw. This test checks the
    // throw. Once we figure out a way to work with underdetermined systems,
    // the test will have to be changed.
    TEST_THROW( ( DataTransferKit::MovingLeastSquaresOperator<
                    DeviceType, RadialBasisFunction, PolynomialBasis>(
                    comm, source_points, target_points, plist ) ),
                DataTransferKit::DataTransferKitException );

    (void)source_points;
    (void)target_points;
    (void)source_values;
    (void)target_values;
    (void)target_values_ref;
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MovingLeastSquaresOperator,
                                   unique_source_point, DeviceType )
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

    TEST_THROW( DataTransferKit::MovingLeastSquaresOperator<DeviceType>(
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

    DataTransferKit::MovingLeastSquaresOperator<DeviceType> mlsop(
        comm, source_points, target_points );

    Kokkos::View<double *, DeviceType> source_values( "in" );
    Kokkos::View<double *, DeviceType> target_values( "out" );

    // violate pre condition of apply
    TEST_THROW( mlsop.apply( source_values, target_values ),
                DataTransferKit::DataTransferKitException );

    Kokkos::realloc( target_values, target_points.extent( 0 ) );
    Kokkos::realloc( source_values, source_points.extent( 0 ) );
    if ( comm_rank == 0 )
    {
        auto source_values_host = Kokkos::create_mirror_view( source_values );
        source_values_host( 0 ) = 255.;
        Kokkos::deep_copy( source_values, source_values_host );
    }

    mlsop.apply( source_values, target_values );

    auto target_values_host = Kokkos::create_mirror_view( target_values );
    Kokkos::deep_copy( target_values_host, target_values );
    std::vector<double> target_values_ref = {255.};
    TEST_COMPARE_ARRAYS( target_values_host, target_values_ref );
}

// Include the test macros.
#include "DataTransferKitMeshfree_ETIHelperMacros.h"

using Wendland0 =
    DataTransferKit::RadialBasisFunction<DataTransferKit::Wendland<0>>;
using Wendland2 =
    DataTransferKit::RadialBasisFunction<DataTransferKit::Wendland<2>>;
using Wendland6 =
    DataTransferKit::RadialBasisFunction<DataTransferKit::Wendland<6>>;
using Constant3 =
    DataTransferKit::MultivariatePolynomialBasis<DataTransferKit::Constant, 3>;
using Linear3 =
    DataTransferKit::MultivariatePolynomialBasis<DataTransferKit::Linear, 3>;
using Quadratic3 =
    DataTransferKit::MultivariatePolynomialBasis<DataTransferKit::Quadratic, 3>;

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    using DeviceType##NODE = typename NODE::device_type;                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(                                      \
        MovingLeastSquaresOperator, same_npoints_and_basis, DeviceType##NODE,  \
        Wendland0, Constant3 )                                                 \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(                                      \
        MovingLeastSquaresOperator, same_npoints_and_basis, DeviceType##NODE,  \
        Wendland0, Linear3 )                                                   \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(                                      \
        MovingLeastSquaresOperator, same_npoints_and_basis, DeviceType##NODE,  \
        Wendland0, Quadratic3 )                                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MovingLeastSquaresOperator,          \
                                          constant_basis, DeviceType##NODE,    \
                                          Wendland0 )                          \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MovingLeastSquaresOperator,          \
                                          constant_basis, DeviceType##NODE,    \
                                          Wendland2 )                          \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MovingLeastSquaresOperator,          \
                                          constant_basis, DeviceType##NODE,    \
                                          Wendland6 )                          \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MovingLeastSquaresOperator,          \
                                          underdetermined, DeviceType##NODE,   \
                                          Wendland0, Linear3 )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
