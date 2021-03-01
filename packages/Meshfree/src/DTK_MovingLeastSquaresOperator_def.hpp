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

#ifndef DTK_MOVING_LEAST_SQUARES_OPERATOR_DEF_HPP
#define DTK_MOVING_LEAST_SQUARES_OPERATOR_DEF_HPP

#include <ArborX.hpp>
#include <DTK_DBC.hpp>
#include <DTK_DetailsMovingLeastSquaresOperatorImpl.hpp>
#include <DTK_DetailsNearestNeighborOperatorImpl.hpp> // fetch
#include <DTK_DetailsUtils.hpp>

namespace DataTransferKit
{

template <typename DeviceType, typename CompactlySupportedRadialBasisFunction,
          typename PolynomialBasis>
MovingLeastSquaresOperator<DeviceType, CompactlySupportedRadialBasisFunction,
                           PolynomialBasis>::
    MovingLeastSquaresOperator(
        MPI_Comm comm,
        Kokkos::View<Coordinate const **, DeviceType> source_points,
        Kokkos::View<Coordinate const **, DeviceType> target_points )
    : _comm( comm )
    , _n_source_points( source_points.extent( 0 ) )
    , _offset( "offset", 0 )
    , _ranks( "ranks", 0 )
    , _indices( "indices", 0 )
    , _coeffs( "polynomial_coefficients", 0 )
{
    DTK_REQUIRE( source_points.extent_int( 1 ) ==
                 target_points.extent_int( 1 ) );
    // FIXME for now let's assume 3D
    DTK_REQUIRE( source_points.extent_int( 1 ) == 3 );

    using ExecutionSpace = typename DeviceType::execution_space;
    using MemorySpace = typename DeviceType::memory_space;
    // Build distributed search tree over the source points.
    ArborX::DistributedTree<MemorySpace> search_tree( _comm, ExecutionSpace{},
                                                      source_points );
    DTK_CHECK( !search_tree.empty() );

    // For each target point, query the n_neighbors points closest to the
    // target.
    auto queries =
        Details::MovingLeastSquaresOperatorImpl<DeviceType>::makeKNNQueries(
            target_points, PolynomialBasis::size );

    // Perform the actual search.
    using PairIndexRank = Kokkos::pair<int, int>;
    Kokkos::View<PairIndexRank *, DeviceType> index_rank( "index_rank", 0 );
    search_tree.query( ExecutionSpace{}, queries, index_rank, _offset );
    // Split the pair
    Details::splitIndexRank( index_rank, _indices, _ranks );

    // Retrieve the coordinates of all source points that met the predicates.
    // NOTE: This is the last collective.
    source_points = Details::NearestNeighborOperatorImpl<DeviceType>::fetch(
        _comm, _ranks, _indices, source_points );

    // Transform source points
    source_points = Details::MovingLeastSquaresOperatorImpl<
        DeviceType>::transformSourceCoordinates( source_points, _offset,
                                                 target_points );
    target_points = Kokkos::View<Coordinate **, DeviceType>( "empty", 0, 0 );

    // Build P (vandermonde matrix)
    // P is a single 1D storage for multiple P_i matrices. Each matrix is of
    // size (#source_points_for_specific_target_point, basis_size)
    auto p =
        Details::MovingLeastSquaresOperatorImpl<DeviceType>::computeVandermonde(
            source_points, PolynomialBasis() );

    // To build the radial basis function, we need to define the radius of the
    // radial basis function. Since we use kNN, we need to compute the radius.
    // We only need the coordinates of the source points because of the
    // transformation of the coordinates.
    auto radius =
        Details::MovingLeastSquaresOperatorImpl<DeviceType>::computeRadius(
            source_points, _offset );

    // Build phi (weight matrix)
    auto phi =
        Details::MovingLeastSquaresOperatorImpl<DeviceType>::computeWeights(
            source_points, radius, CompactlySupportedRadialBasisFunction() );

    // Build A (moment matrix)
    auto a =
        Details::MovingLeastSquaresOperatorImpl<DeviceType>::computeMoments(
            _offset, p, phi );

    // TODO: it is computationally unnecessary to compute the pseudo-inverse as
    // MxM (U*E^+*V) as it will later be just used to do MxV. We could instead
    // return the (U,E^+,V) and do the MxV multiplication. But for now, it's OK.
    auto t = Details::MovingLeastSquaresOperatorImpl<DeviceType>::invertMoments(
        a, PolynomialBasis::size );
    auto inv_a = std::get<0>( t );

    // std::get<1>(t) returns the number of undetermined system. However, this
    // is not enough to know if we will lose order of accuracy. For example, if
    // all the points are aligned, the system will be underdetermined. However
    // this is not a problem if we found at least three points since this is
    // enough to define a quadratic function. Therefore, not only we need to
    // know the rank deficiency but also the dimension of the problem.

    // NOTE: This assumes that the polynomial basis evaluated at {0,0,0} is
    // going to be [1, 0, 0, ..., 0]^T.
    _coeffs = Details::MovingLeastSquaresOperatorImpl<
        DeviceType>::computePolynomialCoefficients( _offset, inv_a, p, phi,
                                                    PolynomialBasis::size );
}

template <typename DeviceType, typename CompactlySupportedRadialBasisFunction,
          typename PolynomialBasis>
void MovingLeastSquaresOperator<
    DeviceType, CompactlySupportedRadialBasisFunction, PolynomialBasis>::
    apply( Kokkos::View<double const *, DeviceType> source_values,
           Kokkos::View<double *, DeviceType> target_values ) const
{
    // Precondition: check that the source and the target are properly sized
    DTK_REQUIRE( source_values.extent( 0 ) == _n_source_points );
    DTK_REQUIRE( target_values.extent( 0 ) == _offset.extent( 0 ) - 1 );

    // Retrieve values for all source points
    source_values = Details::NearestNeighborOperatorImpl<DeviceType>::fetch(
        _comm, _ranks, _indices, source_values );

    // Apply A-1 (P^T phi)
    auto new_target_values = Details::MovingLeastSquaresOperatorImpl<
        DeviceType>::computeTargetValues( _offset, _coeffs, source_values );

    Kokkos::deep_copy( target_values, new_target_values );
}

} // end namespace DataTransferKit

// Explicit instantiation macro
#define DTK_MOVING_LEAST_SQUARES_OPERATOR_INSTANT( NODE )                      \
    template class MovingLeastSquaresOperator<typename NODE::device_type>;     \
    template class MovingLeastSquaresOperator<                                 \
        typename NODE::device_type, Wendland<0>,                               \
        MultivariatePolynomialBasis<Quadratic, 3>>;

#endif
