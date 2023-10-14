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
#include "DTK_Core.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
namespace
{ // anonymous

// Whether one of the Tpetra::initialize() functions has been called before.
bool dtkIsInitialized = false;

// Whether DTK initialized Kokkos. DTK's finalize() only finalizes
// Kokkos if it initialized Kokkos. Otherwise, something else
// initialized Kokkos and is responsible for finalizing it.
bool dtkInitializedKokkos = false;

// Initialize Kokkos, if it needs initialization.
template <typename... Args>
void initKokkos( Args &&... args )
{
    if ( !dtkInitializedKokkos )
    {
        // Kokkos doesn't have a global is_initialized().  However,
        // Kokkos::initialize() always initializes the default execution
        // space, so it suffices to check whether that was initialized.
        const bool kokkosIsInitialized = Kokkos::is_initialized();

        if ( !kokkosIsInitialized )
        {
            // Kokkos will remove all arguments Kokkos recognizes which start
            // with '--kokkos' (e.g.,--kokkos-threads)
            Kokkos::initialize( std::forward<Args>( args )... );
            dtkInitializedKokkos = true;
        }
    }

    const bool kokkosIsInitialized = Kokkos::is_initialized();

    if ( !kokkosIsInitialized )
        throw DataTransferKitException( "At the end of initKokkos, Kokkos"
                                        " is not initialized. Please report"
                                        " this bug to the DTK developers." );
}

} // namespace

template <typename... Args>
void initialize( Args &&... args )
{
    if ( !dtkIsInitialized )
        initKokkos( std::forward<Args>( args )... );
    dtkIsInitialized = true;
}

bool isInitialized() { return dtkIsInitialized; }

void finalize()
{
    if ( !dtkIsInitialized )
        return;

    // DTK should only finalize Kokkos if it initialized it
    if ( dtkInitializedKokkos )
        Kokkos::finalize();

    dtkIsInitialized = false;
}

// ETI for initialize
template void initialize<int &, char **&>( int &argc, char **&argv );
template <>
void initialize<int *, char ***>( int *&&argc, char ***&&argv )
{
    initialize( *argc, *argv );
}
template void initialize<>();
} // namespace DataTransferKit
