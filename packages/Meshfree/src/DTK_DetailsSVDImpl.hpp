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

#ifndef DTK_DETAILS_SVD_IMPL_HPP
#define DTK_DETAILS_SVD_IMPL_HPP

#include <DTK_KokkosHelpers.hpp>

#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_Core.hpp>

#include <cmath>

namespace DataTransferKit
{
namespace Details
{

// The original version of this functor was taken from Trilinos mini-tensor
// package. It was adapted to work in a batched mode where matrices are given
// in a flat 1D array. It also explicitly solves 2x2 singular-value
// decomposition (svd) problems.
template <typename DeviceType>
struct SVDFunctor
{
  public:
    using ExecutionSpace = typename DeviceType::execution_space;

    using matrices_type = Kokkos::View<double *, DeviceType>;
    using shared_matrix =
        Kokkos::View<double **, typename ExecutionSpace::scratch_memory_space,
                     Kokkos::MemoryUnmanaged>;
    using matrix_2x2_type = Kokkos::Array<Kokkos::Array<double, 2>, 2>;

  public:
    SVDFunctor( int n, typename matrices_type::const_type As,
                matrices_type pseudoAs )
        : _n( n )
        , _As( As )
        , _pseudoAs( pseudoAs )
    {
    }

    KOKKOS_INLINE_FUNCTION
    void givens_left( shared_matrix &A, double c, double s, int i, int k ) const
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
    void givens_right( shared_matrix &A, double c, double s, int i,
                       int k ) const
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
    void trans_2x2( const matrix_2x2_type A, matrix_2x2_type &B ) const
    {
        B = {{{{A[0][0], A[1][0]}}, {{A[0][1], A[1][1]}}}};
    }

    KOKKOS_INLINE_FUNCTION
    void trans_nxn( const shared_matrix A, shared_matrix &B ) const
    {
        auto n = A.extent( 0 );

        for ( int i = 0; i < n; i++ )
            for ( int j = 0; j < n; j++ )
                B( i, j ) = A( j, i );
    }

    KOKKOS_INLINE_FUNCTION
    void mult_2x2( const matrix_2x2_type A, const matrix_2x2_type B,
                   matrix_2x2_type &C ) const
    {
        C = {{{{A[0][0] * B[0][0] + A[0][1] * B[1][0],
                A[0][0] * B[0][1] + A[0][1] * B[1][1]}},
              {{A[1][0] * B[0][0] + A[1][1] * B[1][0],
                A[1][0] * B[0][1] + A[1][1] * B[1][1]}}}};
    }

    KOKKOS_INLINE_FUNCTION
    void svd_2x2( const matrix_2x2_type A, matrix_2x2_type &U,
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
    void argmax_off_diagonal( typename shared_matrix::const_type A, int &p,
                              int &q ) const
    {
        const auto n = A.extent_int( 0 );

        p = q = -1;
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
    double norm_F_wo_diag( typename shared_matrix::const_type A ) const
    {
        const auto n = A.extent_int( 0 );

        double norm = 0.0;
        for ( int i = 0; i < n; i++ )
            for ( int j = 0; j < n; j++ )
                norm += ( ( i != j ) ? A( i, j ) * A( i, j ) : 0 );

        return std::sqrt( norm );
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type &thread,
        size_t &num_underdetermined ) const
    {
        // FIXME: This only works when team size is 1. Otherwise, we either have
        // to do
        //   matrix_id = thread.league_rank()*thread.team_size() +
        //   thread.team_rank()
        // or do some shared work for a team.
        int matrix_id = thread.league_rank();

        // TODO: This code (for getting A and pseudoA) can be updated later to
        // work with offsets so that we can solve for matrices of different
        // sizes. However, it is unclear what the best batched approach is.
        // It could be that instead the matrices should be pre-sorted by size.
        auto A = Kokkos::subview(
            _As, Kokkos::make_pair(
                     static_cast<size_t>( matrix_id * _n * _n ),
                     static_cast<size_t>( ( matrix_id + 1 ) * _n * _n ) ) );
        auto pseudoA = Kokkos::subview(
            _pseudoAs,
            Kokkos::make_pair(
                static_cast<size_t>( matrix_id * _n * _n ),
                static_cast<size_t>( ( matrix_id + 1 ) * _n * _n ) ) );

        // Allocate (from scratch) and initialize
        shared_matrix E( thread.team_shmem(), _n, _n );
        shared_matrix U( thread.team_shmem(), _n, _n );
        shared_matrix V( thread.team_shmem(), _n, _n );

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
        auto tol = Kokkos::ArithTraits<double>::epsilon();

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
            auto sr = ( KokkosHelpers::sgn( R[0][1] ) ==
                        KokkosHelpers::sgn( R[1][0] ) )
                          ? -R[0][1]
                          : R[0][1];

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
        num_underdetermined += local_undetermined;
    }

    // amount of shared memory
    size_t team_shmem_size( int team_size ) const
    {
        // If each thread of a team gets its own matrix, we multiply by
        // team_size
        return team_size * 3 * shared_matrix::shmem_size( _n, _n ); // U, E, V
    }

  private:
    int _n;
    typename matrices_type::const_type _As;
    matrices_type _pseudoAs;
};

} // end namespace Details
} // end namespace DataTransferKit

#endif
