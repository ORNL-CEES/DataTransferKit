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

#ifndef DTK_POINTCLOUDPROBLEMGENERATOR_HPP
#define DTK_POINTCLOUDPROBLEMGENERATOR_HPP

#include "DTK_ConfigDefs.hpp"
#include "DTK_Types.h"

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
template <class Scalar, class SourceDevice, class TargetDevice>
class PointCloudProblemGenerator
{
  public:
    virtual ~PointCloudProblemGenerator() = default;

    /*!
     * \brief Create a problem where all points are uniquely owned (i.e. no
     * ghosting)
     *
     * \param src_coords Coordinates of the source points. Layout:
     * (point,dim).
     *
     * \param src_field Multi-component field defined on the source points. At
     * a minimum will be allocated and filled with zeros. Some implementations
     * may choose to fill this view with non-zero data. Layout: (point,comp).
     *
     * \param tgt_coords Coordinates of the target points. Layout:
     * (point,dim).
     *
     * \param tgt_field Multi-component field defined on the target points. At
     * a minimum will be allocated and filled with zeros. Some implementations
     * may choose to fill this view with non-zero data corresponding to the
     * expected result of transferring the src_field from the source points to
     * the target points. Layout: (point,comp).
     */
    virtual void createUniquelyOwnedProblem(
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, SourceDevice>
            &src_coords,
        Kokkos::View<Scalar **, Kokkos::LayoutLeft, SourceDevice> &src_field,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, TargetDevice>
            &tgt_coords,
        Kokkos::View<Scalar **, Kokkos::LayoutLeft, TargetDevice>
            &tgt_field ) = 0;

    /*!
     * \brief Create a general problem where points may exist on multiple
     * processors. Points have a unique global id.
     *
     * \param src_coords Coordinates of the source points. Layout:
     * (point,dim).
     *
     * \param src_gids Global ids of the source points. A global id will
     * appear only once on a given processor but may appear on multiple
     * processors. Layout: (point).
     *
     * \param src_field Multi-component field defined on the source points. At
     * a minimum will be allocated and filled with zeros. Some implementations
     * may choose to fill this view with non-zero data. Layout: (point,comp).
     *
     * \param tgt_coords Coordinates of the target points. Layout:
     * (point,dim).
     *
     * \param tgt_gids Global ids of the target points. A global id will
     * appear only once on a given processor but may appear on multiple
     * processors. Layout: (point).
     *
     * \param tgt_field Multi-component field defined on the target points. At
     * a minimum will be allocated and filled with zeros. Some implementations
     * may choose to fill this view with non-zero data corresponding to the
     * expected result of transferring the src_field from the source points to
     * the target points. Layout: (point,comp).
     */
    virtual void createGhostedProblem(
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, SourceDevice>
            &src_coords,
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, SourceDevice>
            &src_gids,
        Kokkos::View<Scalar **, Kokkos::LayoutLeft, SourceDevice> &src_field,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, TargetDevice>
            &tgt_coords,
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, TargetDevice>
            &tgt_gids,
        Kokkos::View<Scalar **, Kokkos::LayoutLeft, TargetDevice>
            &tgt_field ) = 0;
};

//---------------------------------------------------------------------------//

} // namespace DataTransferKit

#endif // end  DTK_POINTCLOUDPROBLEMGENERATOR_HPP
