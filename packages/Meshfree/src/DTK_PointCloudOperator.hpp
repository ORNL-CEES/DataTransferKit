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

#ifndef DTK_POINT_CLOUD_OPERATOR_DECL_HPP
#define DTK_POINT_CLOUD_OPERATOR_DECL_HPP

#include <DTK_ConfigDefs.hpp>

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{
/**
 * Base class for the MeshFree methods.
 */
template <typename DeviceType>
class PointCloudOperator
{
  public:
    virtual ~PointCloudOperator() = default;

    /**
     * Compute the values of a field at the target points given the values at
     * the source points.
     */
    virtual void
    apply( Kokkos::View<double const *, DeviceType> source_values,
           Kokkos::View<double *, DeviceType> target_values ) const = 0;
};

} // end namespace DataTransferKit

#endif
