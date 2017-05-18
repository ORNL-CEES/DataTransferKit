/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#ifndef DTK_COARSEGLOBALSEARCH_DEF_HPP
#define DTK_COARSEGLOBALSEARCH_DEF_HPP

#include <Teuchos_CommHelpers.hpp>

#include <Tpetra_Distributor.hpp>

#include "DTK_CoarseGlobalSearch.hpp"
#include "DTK_DetailsAlgorithms.hpp"

namespace DataTransferKit
{
template <typename NO>
CoarseGlobalSearch<NO>::CoarseGlobalSearch(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
    Kokkos::View<Box const *, DeviceType> local_boxes )
    : _comm( comm )
{
    using ExecutionSpace = typename DeviceType::execution_space;

    int comm_rank = _comm->getRank();
    int comm_size = _comm->getSize();

    _subdomain_boxes =
        Kokkos::View<Box *, DeviceType>( "subdomain_boxes", comm_size );

    // Assemble the local domain bounding box
    int const n = local_boxes.extent( 0 );
    Details::ExpandBoxWithBoxFunctor<DeviceType> functor( local_boxes );
    Kokkos::parallel_reduce( "calculate_subdomain_bounding",
                             Kokkos::RangePolicy<ExecutionSpace>( 0, n ),
                             functor, _subdomain_boxes( comm_rank ) );

    auto host_subdomain_boxes = Kokkos::create_mirror_view( _subdomain_boxes );
    Kokkos::deep_copy( host_subdomain_boxes, _subdomain_boxes );

    // Gather the bounding boxes from all domains.
    // FIXME: can we do this in place?
    Teuchos::Array<double> all_bounds( 6 * comm_size );
    Teuchos::gatherAll<int, double>(
        *_comm, 6,
        reinterpret_cast<double *>( &host_subdomain_boxes( comm_rank ) ),
        comm_size * 6, all_bounds.getRawPtr() );

    // Extract the bounding boxes.
    for ( int i = 0; i < comm_size; ++i )
    {
        int id = 6 * i;
        host_subdomain_boxes( i ) = {all_bounds[id],     all_bounds[id + 1],
                                     all_bounds[id + 2], all_bounds[id + 3],
                                     all_bounds[id + 4], all_bounds[id + 5]};
    }

    Kokkos::deep_copy( _subdomain_boxes, host_subdomain_boxes );
}
}

// Explicit instantiation macro
#define DTK_COARSEGLOBALSEARCH_INSTANT( NODE )                                 \
    template class CoarseGlobalSearch<NODE>;

#endif
