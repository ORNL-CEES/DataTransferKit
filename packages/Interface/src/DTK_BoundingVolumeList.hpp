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
/*!
 * \file DTK_BoundingVolumeList.hpp
 * \brief BoundingVolume list.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BOUNDINGVOLUMELIST_HPP
#define DTK_BOUNDINGVOLUMELIST_HPP

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class BoundingVolumeList.
 *
 * \brief Trivially-copyable bounding volume list.
 *
 * \tparam ViewProperties Properties of the contained Kokkos views.
 */
template <class... ViewProperties>
class BoundingVolumeList
{
  public:
    //! View traits.
    using ViewTraits = typename Kokkos::ViewTraits<int, ViewProperties...>;

    //! The coordinates of the bounding volumes that are locally-owned by this
    //! MPI rank. This view is rank-3 and should be sized as (number of
    //! bounding volumes, spatial dimension, 2 [for low and high coordinate]).
    Kokkos::View<Coordinate * * [2], ViewProperties...> bounding_volumes;
};

//---------------------------------------------------------------------------//

} // namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_BOUNDINGVOLUMELIST_HPP

//---------------------------------------------------------------------------//
// end DTK_BoundingVolumeList.hpp
//---------------------------------------------------------------------------//
