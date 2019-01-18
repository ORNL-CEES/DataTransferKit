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

#ifndef DTK_NEAREST_NEIGHBOR_OPERATOR_DEF_HPP
#define DTK_NEAREST_NEIGHBOR_OPERATOR_DEF_HPP

#include <DTK_DetailsNearestNeighborOperatorImpl.hpp>
#include <DTK_DistributedSearchTree.hpp>

namespace DataTransferKit
{

template <typename DeviceType>
NearestNeighborOperator<DeviceType>::NearestNeighborOperator(
    MPI_Comm comm, Kokkos::View<Coordinate const **, DeviceType> source_points,
    Kokkos::View<Coordinate const **, DeviceType> target_points )
    : _comm( comm )
    , _indices( "indices" )
    , _ranks( "ranks" )
    , _size( source_points.extent_int( 0 ) )
{
    // NOTE: instead of checking the pre-condition that there is at least one
    // source point passed to one of the rank, we let the tree handle the
    // communication and just check that the tree is not empty.

    // Build distributed search tree over the source points.
    DistributedSearchTree<DeviceType> search_tree( _comm, source_points );

    // Tree must have at least one leaf, otherwise it makes little sense to
    // perform the search for nearest neighbors.
    DTK_CHECK( !search_tree.empty() );

    // Query nearest neighbor for all target points.
    auto nearest_queries = Details::NearestNeighborOperatorImpl<
        DeviceType>::makeNearestNeighborQueries( target_points );

    // Perform the actual search.
    Kokkos::View<int *, DeviceType> indices( "indices" );
    Kokkos::View<int *, DeviceType> offset( "offset" );
    Kokkos::View<int *, DeviceType> ranks( "ranks" );
    search_tree.query( nearest_queries, indices, offset, ranks );

    // Check post-condition that we did find a nearest neighbor to all target
    // points.
    DTK_ENSURE( lastElement( offset ) == target_points.extent_int( 0 ) );

    // Save results.
    // NOTE: we don't bother keeping `offset` around since it is just `[0, 1, 2,
    // ..., n_target_poins]`
    _indices = indices;
    _ranks = ranks;
}

template <typename DeviceType>
void NearestNeighborOperator<DeviceType>::apply(
    Kokkos::View<double const *, DeviceType> source_values,
    Kokkos::View<double *, DeviceType> target_values ) const
{
    // Precondition: check that the source and target are properly sized
    DTK_REQUIRE( _indices.extent( 0 ) == target_values.extent( 0 ) );
    DTK_REQUIRE( _size == source_values.extent_int( 0 ) );

    auto values = Details::NearestNeighborOperatorImpl<DeviceType>::fetch(
        _comm, _ranks, _indices, source_values );

    Kokkos::deep_copy( target_values, values );
}

} // namespace DataTransferKit

// Explicit instantiation macro
#define DTK_NEARESTNEIGHBOROPERATOR_INSTANT( NODE )                            \
    template class NearestNeighborOperator<typename NODE::device_type>;

#endif
