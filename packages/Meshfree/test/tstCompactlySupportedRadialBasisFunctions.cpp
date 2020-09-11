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

#include <Kokkos_Core.hpp>

#include <DTK_CompactlySupportedRadialBasisFunctions.hpp>

#include <boost/math/tools/polynomial.hpp>
#include <boost/math/tools/rational.hpp>

template <typename DeviceType, typename RadialBasisFunction>
void check_polynomial( boost::math::tools::polynomial<double> const &poly,
                       std::vector<double> const &radii,
                       RadialBasisFunction const &rbf,
                       Teuchos::FancyOStream &out, bool &success )
{
    using ExecutionSpace = typename DeviceType::execution_space;
    int const n = radii.size();
    Kokkos::View<double *, DeviceType> r( "radii", n );
    auto r_host = Kokkos::create_mirror_view( r );
    for ( int i = 0; i < n; ++i )
        r_host( i ) = radii[i];
    Kokkos::deep_copy( r, r_host );

    std::vector<double> p = poly.data();
    std::vector<double> values;
    for ( auto const &x : radii )
        values.push_back(
            boost::math::tools::evaluate_polynomial( p.data(), x, p.size() ) );

    Kokkos::View<double *, DeviceType> v( "values", n );
    Kokkos::parallel_for( "evaluate",
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n ),
                          KOKKOS_LAMBDA( int i ) { v( i ) = rbf( r( i ) ); } );
    Kokkos::fence();
    auto v_host = Kokkos::create_mirror_view( v );
    Kokkos::deep_copy( v_host, v );
    double const relative_tolerance = 1.0e-8;
    TEST_COMPARE_FLOATING_ARRAYS( v_host, values, relative_tolerance );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CompactlySupportedRadialBasisFunctions,
                                   polynomial_rbf, DeviceType )
{
    int const n = 10;
    std::vector<double> r( n );
    for ( int i = 0; i < n; ++i )
        r[i] = static_cast<double>( i ) / n;

    check_polynomial<DeviceType>(
        boost::math::tools::pow(
            boost::math::tools::polynomial<double>{1.0, -1.0}, 2 ),
        r, DataTransferKit::Wendland<0>(), out, success );

    check_polynomial<DeviceType>(
        boost::math::tools::pow(
            boost::math::tools::polynomial<double>{1.0, -1.0}, 4 ) *
            boost::math::tools::polynomial<double>{1.0, 4.0},
        r, DataTransferKit::Wendland<2>(), out, success );

    check_polynomial<DeviceType>(
        boost::math::tools::pow(
            boost::math::tools::polynomial<double>{1.0, -1.0}, 6 ) *
            boost::math::tools::polynomial<double>{3.0, 18.0, 35.0},
        r, DataTransferKit::Wendland<4>(), out, success );

    check_polynomial<DeviceType>(
        boost::math::tools::pow(
            boost::math::tools::polynomial<double>{1.0, -1.0}, 8 ) *
            boost::math::tools::polynomial<double>{1.0, 8.0, 25.0, 32.0},
        r, DataTransferKit::Wendland<6>(), out, success );

    check_polynomial<DeviceType>(
        boost::math::tools::pow(
            boost::math::tools::polynomial<double>{1.0, -1.0}, 4 ) *
            boost::math::tools::polynomial<double>{4.0, 16.0, 12.0, 3.0},
        r, DataTransferKit::Wu<2>(), out, success );

    check_polynomial<DeviceType>(
        boost::math::tools::pow(
            boost::math::tools::polynomial<double>{1.0, -1.0}, 6 ) *
            boost::math::tools::polynomial<double>{6.0, 36.0, 82.0, 72.0, 30.0,
                                                   5.0},
        r, DataTransferKit::Wu<4>(), out, success );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CompactlySupportedRadialBasisFunctions,
                                   wrap_rbf, DeviceType )
{
    struct X
    {
        KOKKOS_INLINE_FUNCTION double operator()( double x ) const { return x; }
    };
    DataTransferKit::RadialBasisFunction<X> rbf( 2. );
    TEST_EQUALITY( rbf( 1. ), .5 );
    TEST_EQUALITY( rbf( 2. ), 1. );
    TEST_EQUALITY( rbf( 4. ), 2. );
}

// Include the test macros.
#include "DataTransferKit_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    using DeviceType##NODE = typename NODE::device_type;                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                      \
        CompactlySupportedRadialBasisFunctions, polynomial_rbf,                \
        DeviceType##NODE )                                                     \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                      \
        CompactlySupportedRadialBasisFunctions, wrap_rbf, DeviceType##NODE )
// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
