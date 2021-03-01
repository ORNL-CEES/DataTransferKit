/****************************************************************************
 * Copyright (c) 2012-2021 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef DTK_DETAILS_UTILS_HPP
#define DTK_DETAILS_UTILS_HPP

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
namespace Details
{
template <typename DeviceType>
void splitIndexRank(
    Kokkos::View<Kokkos::pair<int, int> *, DeviceType> index_rank,
    Kokkos::View<int *, DeviceType> &indices,
    Kokkos::View<int *, DeviceType> &ranks )
{
    using ExecutionSpace = typename DeviceType::execution_space;

    auto const n_pairs = index_rank.extent( 0 );
    Kokkos::realloc( indices, n_pairs );
    Kokkos::realloc( ranks, n_pairs );
    Kokkos::parallel_for( DTK_MARK_REGION( "split_pairs" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_pairs ),
                          KOKKOS_LAMBDA( int i ) {
                              indices( i ) = index_rank( i ).first;
                              ranks( i ) = index_rank( i ).second;
                          } );
}
} // namespace Details
} // namespace DataTransferKit
#endif
