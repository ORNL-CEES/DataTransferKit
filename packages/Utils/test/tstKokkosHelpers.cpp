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

#include <DTK_KokkosHelpers.hpp> // isFinite, infinity, isNan

#include <Teuchos_UnitTestHarness.hpp>

#include <cfloat> // DBL_MIN
#include <cmath>  // NAN, INFINITY

// is_nan and is_finite unit tests below adapted from examples on
// cppreference.com

TEUCHOS_UNIT_TEST( KokkosHelpers, is_nan )
{
    using DataTransferKit::KokkosHelpers::isNan;
    TEST_ASSERT( isNan( NAN ) );
    TEST_ASSERT( !isNan( INFINITY ) );
    TEST_ASSERT( !isNan( 0.0 ) );
    TEST_ASSERT( !isNan( DBL_MIN / 2.0 ) );
    TEST_ASSERT( isNan( 0.0 / 0.0 ) );
    TEST_ASSERT( isNan( INFINITY - INFINITY ) );
}

TEUCHOS_UNIT_TEST( KokkosHelpers, is_finite )
{
    using DataTransferKit::KokkosHelpers::isFinite;
    TEST_ASSERT( !isFinite( NAN ) );
    TEST_ASSERT( !isFinite( INFINITY ) );
    TEST_ASSERT( isFinite( 0.0 ) );
    TEST_ASSERT( isFinite( DBL_MIN / 2.0 ) );
    TEST_ASSERT( isFinite( 1.0 ) );
    TEST_ASSERT( !isFinite( std::exp( 800 ) ) );
    TEST_ASSERT( !isFinite(
        DataTransferKit::KokkosHelpers::ArithTraits<double>::infinity() ) );
}
