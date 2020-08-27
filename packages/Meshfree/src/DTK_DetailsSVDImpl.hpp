/****************************************************************************
 * Copyright (c) 2012-2019 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef DTK_DETAILS_SVD_IMPL_HPP
#define DTK_DETAILS_SVD_IMPL_HPP

#include <ArborX_DetailsKokkosExt.hpp> // ArithmeticTraits

#include <Kokkos_Core.hpp>

#include <cassert>
#include <cmath>

namespace DataTransferKit
{
namespace Details
{

/**
 * Branchless sign function. Return 1 if @param x is greater than zero, 0 if
 * @param x is zero, and -1 if @param x is less than zero.
 */
template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
KOKKOS_INLINE_FUNCTION int sgn( T x )
{
    return ( x > 0 ) - ( x < 0 );
}

// The original version of this functor was taken from Trilinos mini-tensor
// package. It was adapted to work in a batched mode where matrices are given
// in a flat 1D array. It also explicitly solves 2x2 singular-value
// decomposition (svd) problems.
template <typename DeviceType>
struct SVDFunctor
{
  public:
    using ExecutionSpace = typename DeviceType::execution_space;

    // We use 1D view of the matrices here to make it as generic as possible.
    // This should allow for a certain flexibility later, like using different
    // sized matrices, or using some batching.
    // NOTE pass flat::matrix_type and matrix_type by value (or typename
    // matrix_type::const_type) but pass matrix_2x2_type by reference (or const
    // &)
    using flat_matrix_type = Kokkos::View<double *, DeviceType>;
    using matrix_type = Kokkos::View<double **, DeviceType>;
    using matrix_2x2_type = Kokkos::Array<Kokkos::Array<double, 2>, 2>;

  public:
    SVDFunctor( int n, typename flat_matrix_type::const_type As,
                flat_matrix_type pseudoAs, matrix_type aux )
        : _n( n )
        , _As( As )
        , _pseudoAs( pseudoAs )
        , _aux( aux )
    {
    }

    KOKKOS_INLINE_FUNCTION
    void givens_left( matrix_type A, double c, double s, int i, int k ) const
    {
        auto n = A.extent_int( 0 );

        for ( int j = 0; j < n; j++ )
        {
            auto aij = A( i, j );
            auto akj = A( k, j );
            A( i, j ) = c * aij - s * akj;
            A( k, j ) = s * aij + c * akj;
        }
    }

    KOKKOS_INLINE_FUNCTION
    void givens_right( matrix_type A, double c, double s, int i, int k ) const
    {
        auto n = A.extent_int( 0 );

        for ( int j = 0; j < n; ++j )
        {
            auto aji = A( j, i );
            auto ajk = A( j, k );
            A( j, i ) = c * aji - s * ajk;
            A( j, k ) = s * aji + c * ajk;
        }
    }

    KOKKOS_INLINE_FUNCTION
    void trans_2x2( matrix_2x2_type const &A, matrix_2x2_type &B ) const
    {
        B = {{{{A[0][0], A[1][0]}}, {{A[0][1], A[1][1]}}}};
    }

    KOKKOS_INLINE_FUNCTION
    void trans_nxn( typename matrix_type::const_type A, matrix_type &B ) const
    {
        auto n = A.extent( 0 );

        for ( int i = 0; i < n; i++ )
            for ( int j = 0; j < n; j++ )
                B( i, j ) = A( j, i );
    }

    KOKKOS_INLINE_FUNCTION
    void mult_2x2( matrix_2x2_type const &A, matrix_2x2_type const &B,
                   matrix_2x2_type &C ) const
    {
        C = {{{{A[0][0] * B[0][0] + A[0][1] * B[1][0],
                A[0][0] * B[0][1] + A[0][1] * B[1][1]}},
              {{A[1][0] * B[0][0] + A[1][1] * B[1][0],
                A[1][0] * B[0][1] + A[1][1] * B[1][1]}}}};
    }

    KOKKOS_INLINE_FUNCTION
    void svd_2x2( matrix_2x2_type const &A, matrix_2x2_type &U,
                  matrix_2x2_type &E, matrix_2x2_type &V ) const
    {
        matrix_2x2_type At, AAt, AtA;
        trans_2x2( A, At );
        mult_2x2( A, At, AAt );
        mult_2x2( At, A, AtA );

        // Find U such that U*A*A’*U’ = diag
        auto phi = 0.5 * atan2( AAt[0][1] + AAt[1][0], AAt[0][0] - AAt[1][1] );
        auto cphi = cos( phi );
        auto sphi = sin( phi );

        U = {{{{cphi, -sphi}}, {{sphi, cphi}}}};

        // Find W such that W’*A’*A*W = diag
        auto theta =
            0.5 * atan2( AtA[0][1] + AtA[1][0], AtA[0][0] - AtA[1][1] );
        auto ctheta = cos( theta );
        auto stheta = sin( theta );
        matrix_2x2_type W = {{{{ctheta, -stheta}}, {{stheta, ctheta}}}};

        // Find the singular values from U
        auto sum = AAt[0][0] + AAt[1][1];
        auto dif = sqrt( ( AAt[0][0] - AAt[1][1] ) * ( AAt[0][0] - AAt[1][1] ) +
                         4 * AAt[0][1] * AAt[1][0] );
        E = {{{{sqrt( 0.5 * ( sum + dif ) ), 0.}},
              {{0., sqrt( 0.5 * ( sum - dif ) )}}}};

        // Find the correction matrix for the right side (S = U'*A*W)
        matrix_2x2_type Ut, AW, S;
        mult_2x2( A, W, AW );
        trans_2x2( U, Ut );
        mult_2x2( Ut, AW, S );

        // We need copysign here to work with singular systems. Using the
        // regular sgn will produce a singular C which would lead to singular
        // V.
        matrix_2x2_type C = {{{{std::copysign( 1., S[0][0] ), 0.0}},
                              {{0.0, std::copysign( 1., S[1][1] )}}}};

        mult_2x2( W, C, V );
    }

    KOKKOS_INLINE_FUNCTION
    void argmax_off_diagonal( typename matrix_type::const_type A, int &p,
                              int &q ) const
    {
        const auto n = A.extent_int( 0 );

        p = -1;
        q = -1;
        double max = -1;

        for ( int i = 0; i < n; i++ )
            for ( int j = 0; j < n; j++ )
                if ( i != j && std::abs( A( i, j ) ) > max )
                {
                    p = i;
                    q = j;
                    max = std::abs( A( i, j ) );
                }
    }

    KOKKOS_INLINE_FUNCTION
    double norm_F_wo_diag( typename matrix_type::const_type A ) const
    {
        const auto n = A.extent_int( 0 );

        double norm = 0.0;
        for ( int i = 0; i < n; i++ )
            for ( int j = 0; j < n; j++ )
                norm += ( ( i != j ) ? A( i, j ) * A( i, j ) : 0 );

        return std::sqrt( norm );
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const int matrix_id, size_t &num_underdetermined ) const
    {
        // TODO: This code (for getting A and pseudoA) can be updated later
        // to work with offsets so that we can solve for matrices of
        // different sizes. However, it is unclear what the best batched
        // approach is. It could be that instead the matrices should be
        // pre-sorted by size.
        auto A = Kokkos::subview(
            _As, Kokkos::make_pair( matrix_id * _n * _n,
                                    ( matrix_id + 1 ) * _n * _n ) );
        auto pseudoA = Kokkos::subview(
            _pseudoAs, Kokkos::make_pair( matrix_id * _n * _n,
                                          ( matrix_id + 1 ) * _n * _n ) );

        auto E = Kokkos::subview(
            _aux, Kokkos::ALL(),
            Kokkos::make_pair( 3 * matrix_id * _n, 3 * matrix_id * _n + _n ) );
        auto U =
            Kokkos::subview( _aux, Kokkos::ALL(),
                             Kokkos::make_pair( 3 * matrix_id * _n + _n,
                                                3 * matrix_id * _n + 2 * _n ) );
        auto V =
            Kokkos::subview( _aux, Kokkos::ALL(),
                             Kokkos::make_pair( 3 * matrix_id * _n + 2 * _n,
                                                3 * matrix_id * _n + 3 * _n ) );

        for ( int i = 0; i < _n; i++ )
            for ( int j = 0; j < _n; j++ )
            {
                E( i, j ) = A( i * _n + j );
            }
        for ( int i = 0; i < _n; i++ )
            for ( int j = 0; j < _n; j++ )
            {
                U( i, j ) = ( i == j ? 1.0 : 0.0 );
                V( i, j ) = ( i == j ? 1.0 : 0.0 );
            }

        auto norm = norm_F_wo_diag( E );
        auto tol = KokkosExt::ArithmeticTraits::epsilon<double>::value;

        while ( norm > tol )
        {
            // Find largest off-diagonal entry
            int p, q;
            argmax_off_diagonal( E, p, q );
            assert( p != -1 && q != -1 );
            // TODO: it is unclear whether this permutation is necessary.
            if ( p > q )
            {
                auto t = p;
                p = q;
                q = t;
            }

            // Obtain left and right Givens rotations by using 2x2 SVD
            matrix_2x2_type Apq = {
                {{{E( p, p ), E( p, q )}}, {{E( q, p ), E( q, q )}}}};
            matrix_2x2_type L, D, R;

            svd_2x2( Apq, L, D, R );

            auto cl = L[0][0];
            auto sl = L[0][1];
            auto cr = R[0][0];
            auto sr = ( sgn( R[0][1] ) == sgn( R[1][0] ) ) ? -R[0][1] : R[0][1];

            // Apply both Givens rotations to matrices that are converging to
            // singular values and singular vectors
            givens_left( E, cl, sl, p, q );
            givens_right( E, cr, sr, p, q );

            givens_right( U, cl, sl, p, q );
            givens_left( V, cr, sr, p, q );

            norm = norm_F_wo_diag( E );
        }

        // Compute pseudo-inverse (pseudoA = V pseudoE U^T)
        // NOTE: the V stored above is actually V^T, but we don't explicitly
        // transpose it. Instead, we modify the MxM loop below to do (pseudoA =
        // V^T pseudoE U^T)
        size_t local_undetermined = 0;
        for ( int i = 0; i < _n; i++ )
            for ( int j = 0; j < _n; j++ )
            {
                double value = 0;
                for ( int k = 0; k < _n; k++ )
                {
                    // TODO: We use machine tolerance here to indicate that all
                    // diagonal values less than that are considered to be 0. It
                    // is unclear the numerical implications of such approach
                    // for matrices where singular values are small nonzeros.
                    if ( std::abs( E( k, k ) ) >= tol )
                        value += V( k, i ) * U( j, k ) / E( k, k );
                    else
                        local_undetermined = 1;
                }
                pseudoA( i * _n + j ) = value;
            }
        // TODO: when the kernel is switched to multiple threads per team, this
        // should be fixed. For example, could be an atomic update (as local
        // counts are not shared).
        num_underdetermined += local_undetermined;
    }

  private:
    int _n;
    typename flat_matrix_type::const_type _As;
    flat_matrix_type _pseudoAs;
    matrix_type _aux;
};

} // end namespace Details
} // end namespace DataTransferKit

#endif
