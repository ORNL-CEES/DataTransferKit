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

#include <DTK_DBC.hpp>
#include <DTK_DetailsSVDImpl.hpp>

#include <Kokkos_Core.hpp>

#include <random>
#include <set>

template <typename DeviceType>
void check_result( Kokkos::View<double *, DeviceType> matrices,
                   Kokkos::View<double *, DeviceType> inv_matrices,
                   int const n_matrices, int matrix_size,
                   std::set<int> const &rank_deficiency,
                   Teuchos::FancyOStream &out, bool &success )
{
    auto matrices_host = Kokkos::create_mirror_view( matrices );
    Kokkos::deep_copy( matrices_host, matrices );
    auto inv_matrices_host = Kokkos::create_mirror_view( inv_matrices );
    Kokkos::deep_copy( inv_matrices_host, inv_matrices );

    // Multiply the matrices with their inverse and check that the result is the
    // identity matrix.
    int const offset = matrix_size * matrix_size;
    for ( int m = 0; m < n_matrices; ++m )
    {
        std::vector<std::vector<double>> result(
            matrix_size, std::vector<double>( matrix_size, 0. ) );
        for ( int i = 0; i < matrix_size; ++i )
            for ( int j = 0; j < matrix_size; ++j )
                for ( int k = 0; k < matrix_size; ++k )
                    result[i][j] +=
                        inv_matrices_host[m * offset + i * matrix_size + k] *
                        matrices_host[m * offset + k * matrix_size + j];

        // Check that the result is the identity
        double const relative_tolerance = 1e-12;
        for ( int i = 0; i < matrix_size; ++i )
            for ( int j = 0; j < matrix_size; ++j )
            {
                if ( ( i == j ) && ( rank_deficiency.count( i ) == 0 ) )
                {
                    TEST_FLOATING_EQUALITY( result[i][i], 1.,
                                            relative_tolerance );
                }
                else
                {
                    TEST_FLOATING_EQUALITY( result[i][j] + 1, 1.,
                                            relative_tolerance );
                }
            }
    }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SVD, full_rank, DeviceType )
{
    int const n_matrices = 10;
    int const matrix_size = 32;
    int const size = n_matrices * matrix_size * matrix_size;
    Kokkos::View<double *, DeviceType> matrices( "matrices", size );
    Kokkos::View<double *, DeviceType> inv_matrices( "inv_matrices", size );
    // For magic number 3, see comment in
    // DTK_DetailsMovingLeastSquaresOperatorImpl.hpp
    Kokkos::View<double **, DeviceType> aux( "aux", matrix_size,
                                             3 * n_matrices * matrix_size );

    // Fill the matrices
    auto matrices_host = Kokkos::create_mirror_view( matrices );
    std::default_random_engine random_engine;
    std::uniform_real_distribution<double> distribution( -1000, 1000 );
    for ( int i = 0; i < size; ++i )
        matrices_host( i ) = distribution( random_engine );
    Kokkos::deep_copy( matrices, matrices_host );

    DataTransferKit::Details::SVDFunctor<DeviceType> svd_functor(
        matrix_size, matrices, inv_matrices, aux );
    size_t n_underdetermined = 0;
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::parallel_reduce(
        DTK_MARK_REGION( "compute_svd_inverse" ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, n_matrices ), svd_functor,
        n_underdetermined );

    std::set<int> rank_deficiency;

    check_result( matrices, inv_matrices, n_matrices, matrix_size,
                  rank_deficiency, out, success );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SVD, rank_deficient, DeviceType )
{
    int const n_matrices = 10;
    int const matrix_size = 40;
    int const size = n_matrices * matrix_size * matrix_size;
    Kokkos::View<double *, DeviceType> matrices( "matrices", size );
    Kokkos::View<double *, DeviceType> inv_matrices( "inv_matrices", size );
    // For magic number 3, see comment in
    // DTK_DetailsMovingLeastSquaresOperatorImpl.hpp
    Kokkos::View<double **, DeviceType> aux( "aux", matrix_size,
                                             3 * n_matrices * matrix_size );

    // Fill the matrices
    auto matrices_host = Kokkos::create_mirror_view( matrices );
    std::default_random_engine random_engine;
    std::uniform_real_distribution<double> distribution( -1000, 1000 );
    std::set<int> rank_deficiency;
    rank_deficiency.insert( 3 );
    unsigned int pos = 0;
    for ( int i = 0; i < n_matrices; ++i )
    {
        for ( int j = 0; j < matrix_size; ++j )
        {
            for ( int k = 0; k < matrix_size; ++k )
            {
                if ( rank_deficiency.count( k ) == 0 )
                    matrices_host( pos ) = distribution( random_engine );
                else
                    matrices_host( pos ) = 0.;
                ++pos;
            }
        }
    }
    Kokkos::deep_copy( matrices, matrices_host );

    DataTransferKit::Details::SVDFunctor<DeviceType> svd_functor(
        matrix_size, matrices, inv_matrices, aux );
    size_t n_underdetermined = 0;
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::parallel_reduce(
        DTK_MARK_REGION( "compute_svd_inverse" ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, n_matrices ), svd_functor,
        n_underdetermined );

    TEST_EQUALITY( n_underdetermined, n_matrices );

    check_result( matrices, inv_matrices, n_matrices, matrix_size,
                  rank_deficiency, out, success );
}

// Include the test macros.
#include "DataTransferKit_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    using DeviceType##NODE = typename NODE::device_type;                       \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SVD, full_rank, DeviceType##NODE )   \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( SVD, rank_deficient,                 \
                                          DeviceType##NODE )
// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
