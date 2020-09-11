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

#include <ArborX.hpp>
#include <DTK_MultivariatePolynomialBasis.hpp>

#include <Kokkos_Array.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <vector>

template <typename PolynomialBasis, typename Point>
void checkBasisEvaluation( PolynomialBasis const &b, Point const &x,
                           std::vector<double> const &b_ref, bool &success,
                           Teuchos::FancyOStream &out )
{
    TEST_COMPARE_ARRAYS( b( x ), b_ref );
}

// NOTE: The extra pairs of parentheses below around
// MultivariatePolynomialBasis::size in TEST_EQUALITY() macro are needed.
// Without them, the Teuchos unit testing macro yields the folowing error:
//
//     error: macro "TEST_EQUALITY" passed 3 arguments, but takes just 2

TEUCHOS_UNIT_TEST( MultivariatePolynomialBasis, 2D )
{
    using DataTransferKit::Constant;
    using DataTransferKit::Linear;
    using DataTransferKit::MultivariatePolynomialBasis;
    using DataTransferKit::Quadratic;
    // FIXME
    using Point = Kokkos::Array<double, 2>;

    // (X, Y) -> [ 1 ]
    TEST_EQUALITY( ( MultivariatePolynomialBasis<Constant, 2>::size ), 1 );
    checkBasisEvaluation( MultivariatePolynomialBasis<Constant, 2>(),
                          Point{{0., 0.}}, {1.}, success, out );

    // (X, Y) -> [ 1, X, Y ]
    TEST_EQUALITY( ( MultivariatePolynomialBasis<Linear, 2>::size ), 3 );
    checkBasisEvaluation( MultivariatePolynomialBasis<Linear, 2>(),
                          Point{{0., 1.}}, {{1., 0., 1.}}, success, out );

    // (X, Y) -> [ 1, X, Y, X^2, XY, Y^2 ]
    TEST_EQUALITY( ( MultivariatePolynomialBasis<Quadratic, 2>::size ), 6 );
    checkBasisEvaluation( MultivariatePolynomialBasis<Quadratic, 2>(),
                          Point{{1., 2.}}, {{1., 1., 2., 1., 2., 4.}}, success,
                          out );
}

TEUCHOS_UNIT_TEST( MultivariatePolynomialBasis, 3D )
{
    using ArborX::Point;
    using DataTransferKit::Constant;
    using DataTransferKit::Linear;
    using DataTransferKit::MultivariatePolynomialBasis;
    using DataTransferKit::Quadratic;

    // (X, Y, Z) -> [ 1 ]
    TEST_EQUALITY( ( MultivariatePolynomialBasis<Constant, 3>::size ), 1 );
    checkBasisEvaluation( MultivariatePolynomialBasis<Constant, 3>(),
                          Point{{0., 0., 0.}}, {1.}, success, out );

    // (X, Y, Z) -> [ 1, X, Y, Z ]
    TEST_EQUALITY( ( MultivariatePolynomialBasis<Linear, 3>::size ), 4 );
    checkBasisEvaluation( MultivariatePolynomialBasis<Linear, 3>(),
                          Point{{0., 1., 2.}}, {{1., 0., 1., 2.}}, success,
                          out );

    // (X, Y, Z) -> [ 1, X, Y, Z, X^2, XY, XZ, Y^2, YZ, Z^2 ]
    TEST_EQUALITY( ( MultivariatePolynomialBasis<Quadratic, 3>::size ), 10 );
    checkBasisEvaluation(
        MultivariatePolynomialBasis<Quadratic, 3>(), Point{{0., 1., 2.}},
        {{1., 0., 1., 2., 0., 0., 0., 1., 2., 4.}}, success, out );
}
