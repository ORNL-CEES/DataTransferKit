/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#ifndef DTK_POINTCLOUDPROBLEMGENERATOR_HPP
#define DTK_POINTCLOUDPROBLEMGENERATOR_HPP

#include "DTK_ConfigDefs.hpp"
#include "DTK_Types.h"

#include <Kokkos_View.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
template <class SourceDevice, class TargetDevice>
class PointCloudProblemGenerator
{
  public:
    virtual ~PointCloudProblemGenerator() = default;

    // Create a problem where all points are uniquely owned (i.e. no ghosting)
    virtual void createUniquelyOwnedProblem(
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, SourceDevice>
            &src_coords,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, TargetDevice>
            &tgt_coords ) = 0;

    // Create a general problem where points may exist on multiple
    // processors. Points have a unique global id.
    virtual void createGhostedProblem(
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, SourceDevice>
            &src_coords,
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, SourceDevice>
            &src_gids,
        Kokkos::View<Coordinate **, Kokkos::LayoutLeft, TargetDevice>
            &tgt_coords,
        Kokkos::View<GlobalOrdinal *, Kokkos::LayoutLeft, TargetDevice>
            &tgt_gids ) = 0;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end  DTK_POINTCLOUDPROBLEMGENERATOR_HPP
