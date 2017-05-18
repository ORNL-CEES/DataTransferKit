/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#ifndef DTK_COARSEGLOBALSEARCH_DECL_HPP
#define DTK_COARSEGLOBALSEARCH_DECL_HPP

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

#include "DTK_ConfigDefs.hpp"
#include "DTK_DetailsBox.hpp"

namespace DataTransferKit
{
template <typename NO>
class CoarseGlobalSearch
{
  public:
    using DeviceType = typename NO::device_type;

    // Constructor.
    CoarseGlobalSearch( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                        Kokkos::View<Box const *, DeviceType> local_boxes );

  private:
    // Communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> _comm;

    // Domain bounding boxes.
    Kokkos::View<Box *, DeviceType> _subdomain_boxes;
};
}
#endif
