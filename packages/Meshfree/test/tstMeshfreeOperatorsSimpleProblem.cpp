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

template <typename Operator>
class MeshfreeOperator
{
    using DeviceType = typename Operator::device_type;

  public:
    MeshfreeOperator(
        std::vector<std::array<DataTransferKit::Coordinate, DIM>> const
            &source_points,
        std::vector<std::array<DataTransferKit::Coordinate, DIM>> const
            &target_points );

    void apply( std::vector<double> const &source_values,
                std::vector<double> &target_values );

  private:
    std::shared_ptr<Operator> _op;
};

template <typename Operator>
MeshfreeOperator<Operator>::MeshfreeOperator(
    std::vector<std::array<DataTransferKit::Coordinate, DIM>> const
        &source_points,
    std::vector<std::array<DataTransferKit::Coordinate, DIM>> const
        &target_points )
{
    unsigned int const n_source_points = source_points.size();
    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> sources(
        "source_points", n_source_points, DIM );
    auto sources_host = Kokkos::create_mirror_view( sources );
    for ( unsigned int i = 0; i < n_source_points; ++i )
        for ( unsigned int j = 0; j < DIM; ++j )
            sources_host( i, j ) = source_points[i][j];
    Kokkos::deep_copy( sources, sources_host );

    unsigned int const n_target_points = target_points.size();
    Kokkos::View<DataTransferKit::Coordinate **, DeviceType> targets(
        "target_points", n_target_points, DIM );
    auto targets_host = Kokkos::create_mirror_view( targets );
    for ( unsigned int i = 0; i < n_target_points; ++i )
        for ( unsigned int j = 0; j < DIM; ++j )
            targets_host( i, j ) = target_points[i][j];
    Kokkos::deep_copy( targets, targets_host );

    _op = std::make_shared<Operator>( MPI_COMM_WORLD, sources, targets );
}

template <typename Operator>
void MeshfreeOperator<Operator>::apply(
    std::vector<double> const &source_values,
    std::vector<double> &target_values )
{
    unsigned int const n_source_values = source_values.size();
    Kokkos::View<double *, DeviceType> sources( "source_values",
                                                n_source_values );
    auto sources_host = Kokkos::create_mirror_view( sources );
    for ( unsigned int i = 0; i < n_source_values; ++i )
        sources_host( i ) = source_values[i];
    Kokkos::deep_copy( sources, sources_host );

    unsigned int const n_target_values = target_values.size();
    Kokkos::View<double *, DeviceType> targets( "target_values",
                                                n_target_values );
    auto targets_host = Kokkos::create_mirror_view( targets );
    for ( unsigned int i = 0; i < n_target_values; ++i )
        targets_host( i ) = target_values[i];
    Kokkos::deep_copy( targets, targets_host );

    _op->apply( sources, targets );

    Kokkos::deep_copy( targets_host, targets );
    for ( unsigned int i = 0; i < n_target_values; ++i )
        target_values[i] = targets_host( i );
}

void checkResults( std::vector<double> const &values,
                   std::vector<double> const &references,
                   Teuchos::FancyOStream &out, bool &success )
{
    double const relative_tolerance = 1e-12;
    TEST_COMPARE_FLOATING_ARRAYS( values, references, relative_tolerance );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MeshfreeOperatorSimpleProblem, corner_cases,
                                   Operator )
{
    // single point
    {
        std::vector<std::array<DataTransferKit::Coordinate, DIM>>
            source_points = {{1., 1., 1.}, {2., 2., 2.}};
        std::vector<std::array<DataTransferKit::Coordinate, DIM>>
            target_points = {{1.5, 1.5, 1.5}};
        MeshfreeOperator<Operator> op( source_points, target_points );

        std::vector<double> source_values = {255., 255.};
        std::vector<double> target_values = {0.};
        op.apply( source_values, target_values );

        std::vector<double> ref_values = {255.};
        checkResults( target_values, ref_values, out, success );
    }

    // One source point but no target point.
    {
        std::vector<std::array<DataTransferKit::Coordinate, DIM>>
            source_points = {{1., 1., 1.}};
        std::vector<std::array<DataTransferKit::Coordinate, DIM>>
            target_points = {};
        MeshfreeOperator<Operator> op( source_points, target_points );

        std::vector<double> source_values = {255.};
        std::vector<double> target_values = {};
        op.apply( source_values, target_values );

        std::vector<double> ref_values = {};
        checkResults( target_values, ref_values, out, success );
    }
}

// Include the test macros.
#include "DataTransferKit_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    using DeviceType##NODE = typename NODE::device_type;                       \
    using RadialBasisFunction = DataTransferKit::Wendland<0>;                  \
    using PolynomialBasis =                                                    \
        DataTransferKit::MultivariatePolynomialBasis<DataTransferKit::Linear,  \
                                                     3>;                       \
    using MLS##NODE = DataTransferKit::MovingLeastSquaresOperator<             \
        DeviceType##NODE, RadialBasisFunction, PolynomialBasis>;               \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MeshfreeOperatorSimpleProblem,       \
                                          corner_cases, MLS##NODE )            \
    using Spline##NODE =                                                       \
        DataTransferKit::SplineOperator<DeviceType##NODE, RadialBasisFunction, \
                                        PolynomialBasis>;                      \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MeshfreeOperatorSimpleProblem,       \
                                          corner_cases, Spline##NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
