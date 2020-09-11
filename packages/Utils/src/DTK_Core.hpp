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
 * \file
 * \brief DTK initialization routines.
 */
#ifndef DTK_CORE_HPP
#define DTK_CORE_HPP

#include <Kokkos_Core.hpp>

namespace DataTransferKit
{

// NOTE: trying to follow Tpetra
// (trilinos/packages/tpetra/core/src/Tpetra_Core.hpp)

/*! Initialize DTK
 *
 * Will initialize Kokkos if it was not previously initialized.
 */
template <typename... Args>
void initialize( Args &&... args );

/*! Whether DTK is in initialized state */
bool isInitialized();

/*! Finalize DTK
 *
 * Will finalize Kokkos if it was initialized by DTK.
 */

void finalize();
} // namespace DataTransferKit

#endif // DTK_CORE_HPP
