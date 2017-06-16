/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
/*!
 * \file
 * \brief DTK initialization routines.
 */
#ifndef DTK_CORE_HPP
#define DTK_CORE_HPP

namespace DataTransferKit
{

// NOTE: trying to follow Tpetra
// (trilinos/packages/tpetra/core/src/Tpetra_Core.hpp)

/*! Initialize DTK
 *
 * Will initialize Kokkos if it was not previously initialized.
 */
void initialize( int *argc, char ***argv );

/*! Whether DTK is in initialized state */
bool isInitialized();

/*! Finalize DTK
 *
 * Will finalize Kokkos if it was initialized by DTK.
 */

void finalize();
}

#endif // DTK_CORE_HPP
