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

#ifndef DTK_DETAILS_MOVING_LEAST_SQUARES_OPERATOR_IMPL_HPP
#define DTK_DETAILS_MOVING_LEAST_SQUARES_OPERATOR_IMPL_HPP

#include <DTK_Box.hpp>
#include <DTK_DetailsSVDImpl.hpp>
#include <DTK_DetailsUtils.hpp> // lastElement
#include <DTK_Point.hpp>
#include <DTK_Predicates.hpp>

namespace DataTransferKit
{
namespace Details
{

template <typename DeviceType>
struct MovingLeastSquaresOperatorImpl
{
    using ExecutionSpace = typename DeviceType::execution_space;

    static Kokkos::View<Within *, DeviceType> makeWithinRadiusQueries(
        typename Kokkos::View<Coordinate **, DeviceType>::const_type points,
        double radius )
    {
        auto const n_points = points.extent( 0 );
        Kokkos::View<Within *, DeviceType> queries( "queries", n_points );
        Kokkos::parallel_for(
            DTK_MARK_REGION( "setup_queries" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_points ),
            KOKKOS_LAMBDA( int i ) {
                queries( i ) = within(
                    Point{{points( i, 0 ), points( i, 1 ), points( i, 2 )}},
                    radius );
            } );
        Kokkos::fence();
        return queries;
    }

    static Kokkos::View<double *, DeviceType> computeTargetValues(
        Kokkos::View<int const *, DeviceType> offset,
        Kokkos::View<double const *, DeviceType> polynomial_coeffs,
        Kokkos::View<double const *, DeviceType> source_values

    )
    {
        auto const n_target_points = offset.extent_int( 0 ) - 1;
        Kokkos::View<double *, DeviceType> target_values(
            std::string( "target_" ) + source_values.label(), n_target_points );

        Kokkos::parallel_for(
            DTK_MARK_REGION( "compute_values" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_target_points ),
            KOKKOS_LAMBDA( const int i ) {
                target_values( i ) = 0.;
                for ( int j = offset( i ); j < offset( i + 1 ); ++j )
                    target_values( i ) +=
                        polynomial_coeffs( j ) * source_values( j );
            } );
        Kokkos::fence();

        return target_values;
    }

    static Kokkos::View<Coordinate **, DeviceType> transformSourceCoordinates(
        Kokkos::View<Coordinate const **, DeviceType> source_points,
        Kokkos::View<int const *, DeviceType> offset,
        Kokkos::View<Coordinate const **, DeviceType> target_points )
    {
        auto const n_source_points = source_points.extent( 0 );
        auto const n_target_points = target_points.extent( 0 );

        int const spatial_dim = 3;
        DTK_REQUIRE( source_points.extent_int( 1 ) == spatial_dim );
        DTK_REQUIRE( offset.extent( 0 ) == n_target_points + 1 );

        Kokkos::View<Coordinate **, DeviceType> new_source_points(
            "transformed_source_coords", n_source_points, spatial_dim );
        Kokkos::parallel_for(
            DTK_MARK_REGION( "transform" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_target_points ),
            KOKKOS_LAMBDA( const int i ) {
                for ( int j = offset( i ); j < offset( i + 1 ); j++ )
                    for ( int k = 0; k < spatial_dim; k++ )
                        new_source_points( j, k ) =
                            source_points( j, k ) - target_points( i, k );
            } );

        return new_source_points;
    }

    template <typename RadialBasisFunction>
    static Kokkos::View<double *, DeviceType>
    computeWeights( Kokkos::View<double const **, DeviceType> source_points,
                    RadialBasisFunction const &rbf )
    {
        auto const n_source_points = source_points.extent( 0 );

        int const spatial_dim = 3;
        DTK_REQUIRE( source_points.extent_int( 1 ) == spatial_dim );

        Kokkos::View<double *, DeviceType> phi( "weights", n_source_points );
        Kokkos::parallel_for(
            DTK_MARK_REGION( "compute_weights" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_source_points ),
            KOKKOS_LAMBDA( int i ) {
                phi( i ) = rbf( Details::distance(
                    Point{{source_points( i, 0 ), source_points( i, 1 ),
                           source_points( i, 2 )}},
                    Point{{0., 0., 0.}} ) );
            } );
        Kokkos::fence();
        return phi;
    }

    template <typename PolynomialBasis>
    static Kokkos::View<double *, DeviceType>
    computeVandermonde( Kokkos::View<double const **, DeviceType> points,
                        PolynomialBasis const &polynomial_basis )
    {
        auto const n_points = points.extent( 0 );
        auto constexpr size_polynomial_basis = PolynomialBasis::size();
        Kokkos::View<double *, DeviceType> p(
            "vandermonde", n_points * size_polynomial_basis );
        Kokkos::parallel_for(
            DTK_MARK_REGION( "compute_polynomial_basis" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_points ),
            KOKKOS_LAMBDA( int i ) {
                auto const tmp = polynomial_basis(
                    Point{{points( i, 0 ), points( i, 1 ), points( i, 2 )}} );
                for ( int j = 0; j < size_polynomial_basis; ++j )
                    p( i * size_polynomial_basis + j ) = tmp[j];
            } );
        return p;
    }

    static Kokkos::View<double *, DeviceType>
    computeMoments( Kokkos::View<int const *, DeviceType> offset,
                    Kokkos::View<double const *, DeviceType> p,
                    Kokkos::View<double const *, DeviceType> phi )
    {
        auto const n_target_points = offset.extent_int( 0 ) - 1;
        auto const n_source_points = phi.extent_int( 0 );
        DTK_REQUIRE( n_source_points == lastElement( offset ) );
        // TODO: explain why this is the correct calculation
        auto const size_polynomial_basis =
            n_source_points > 0 ? p.extent_int( 0 ) / n_source_points : 0;
        auto const size_polynomial_basis_squared =
            size_polynomial_basis * size_polynomial_basis;
        Kokkos::View<double *, DeviceType> a(
            "moments", n_target_points * size_polynomial_basis_squared );
        Kokkos::parallel_for(
            DTK_MARK_REGION( "compute_moments" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_target_points ),
            KOKKOS_LAMBDA( int i ) {
                auto p_i = Kokkos::subview(
                    p, Kokkos::make_pair( offset( i ) * size_polynomial_basis,
                                          offset( i + 1 ) *
                                              size_polynomial_basis ) );
                auto phi_i = Kokkos::subview(
                    phi, Kokkos::make_pair( offset( i ), offset( i + 1 ) ) );
                auto a_i = Kokkos::subview(
                    a, Kokkos::make_pair( i * size_polynomial_basis_squared,
                                          ( i + 1 ) *
                                              size_polynomial_basis_squared ) );
                for ( int j = 0; j < size_polynomial_basis; ++j )
                    for ( int k = 0; k < size_polynomial_basis; ++k )
                    {
                        double tmp = 0.;
                        for ( int l = 0; l < offset( i + 1 ) - offset( i );
                              ++l )
                            // Compute value (j,k)
                            tmp += p_i( l * size_polynomial_basis + j ) *
                                   phi_i( l ) *
                                   p_i( l * size_polynomial_basis + k );
                        a_i( j * size_polynomial_basis + k ) = tmp;
                    }
            } );
        Kokkos::fence();

        return a;
    }

    // Matrix pseudo-inversion using SVD
    // Takes in a 1D array of matrices of size NxN, and returns a 1D array of
    // matrices of the same size containing corresponding pseudo-inverses
    static std::tuple<Kokkos::View<double *, DeviceType>, size_t>
    invertMoments( Kokkos::View<double const *, DeviceType> a,
                   const int size_polynomial_basis )
    {
        Kokkos::View<double *, DeviceType> inv_a( "inv_a", a.extent( 0 ) );

        auto num_matrices =
            a.extent( 0 ) / ( size_polynomial_basis * size_polynomial_basis );

        SVDFunctor<DeviceType> svdFunctor( size_polynomial_basis, a, inv_a );
        size_t num_underdetermined = 0;
        Kokkos::parallel_reduce(
            DTK_MARK_REGION( "compute_svd_inverse" ),
            Kokkos::TeamPolicy<ExecutionSpace>( num_matrices, 1 ), svdFunctor,
            num_underdetermined );

        return std::make_tuple( inv_a, num_underdetermined );
    }

    static Kokkos::View<double *, DeviceType> computePolynomialCoefficients(
        Kokkos::View<int const *, DeviceType> offset,
        Kokkos::View<double const *, DeviceType> inv_a,
        Kokkos::View<double const *, DeviceType> p,
        Kokkos::View<double const *, DeviceType> phi,
        const int size_polynomial_basis )
    {
        auto const size_polynomial_basis_squared =
            size_polynomial_basis * size_polynomial_basis;

        auto num_matrices = inv_a.extent( 0 ) / size_polynomial_basis_squared;

        Kokkos::View<double *, DeviceType> coeffs( "polynomial_coeffs",
                                                   phi.extent( 0 ) );

        Kokkos::parallel_for(
            DTK_MARK_REGION( "compute_polynomial_coeffs" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, num_matrices ),
            KOKKOS_LAMBDA( const int i ) {
                auto p_i = Kokkos::subview(
                    p, Kokkos::make_pair( offset( i ) * size_polynomial_basis,
                                          offset( i + 1 ) *
                                              size_polynomial_basis ) );
                auto phi_i = Kokkos::subview(
                    phi, Kokkos::make_pair( offset( i ), offset( i + 1 ) ) );
                auto inv_a_i = Kokkos::subview(
                    inv_a, Kokkos::make_pair(
                               i * size_polynomial_basis_squared,
                               ( i + 1 ) * size_polynomial_basis_squared ) );
                auto coeffs_i = Kokkos::subview(
                    coeffs, Kokkos::make_pair( offset( i ), offset( i + 1 ) ) );

                // coeffs = [1 0 ... 0] * a_inv * p^T * phi
                for ( int k = 0; k < offset( i + 1 ) - offset( i ); k++ )
                {
                    coeffs_i( k ) = 0.;
                    for ( int j = 0; j < size_polynomial_basis; j++ )
                        coeffs_i( k ) +=
                            inv_a_i( 0 * size_polynomial_basis + j ) *
                            p_i( k * size_polynomial_basis + j ) * phi_i( k );
                }
            } );
        return coeffs;
    }
};

} // end namespace Details
} // end namespace DataTransferKit

#endif
