/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#include <DTK_C_API.h>
#include <DTK_C_API_Map.hpp>

#include <cerrno>
#include <set>

//---------------------------------------------------------------------------//
namespace DataTransferKit
{

// We store the reinterpret_cast versions of pointers
static std::set<void *> valid_map_handles;

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

extern "C" {

//---------------------------------------------------------------------------//
DTK_MapHandle DTK_createMap( DTK_ExecutionSpace space, MPI_Comm comm,
                             DTK_UserApplicationHandle source,
                             DTK_UserApplicationHandle target,
                             const char *options )
{
    if ( !DTK_isInitialized() )
    {
        errno = DTK_UNINITIALIZED;
        return nullptr;
    }

    // For demonstration purposes just use the nearest neighbor map.
    auto handle = reinterpret_cast<DTK_MapHandle>(
        DataTransferKit::createMap( space, comm, source, target, options ) );
    DataTransferKit::valid_map_handles.insert( handle );

    errno = DTK_SUCCESS;

    return handle;
}

//---------------------------------------------------------------------------//
bool DTK_isValidMap( DTK_MapHandle handle )
{
    errno = DTK_SUCCESS;
    return DataTransferKit::valid_map_handles.count( handle );
}

//---------------------------------------------------------------------------//
void DTK_applyMap( DTK_MapHandle handle, const char *source_field,
                   const char *target_field )
{
    if ( !DTK_isValidMap( handle ) )
    {
        errno = DTK_INVALID_HANDLE;
        return;
    }

    reinterpret_cast<DataTransferKit::DTK_Map *>( handle )->apply(
        std::string( source_field ), std::string( target_field ) );

    errno = DTK_SUCCESS;
}

//---------------------------------------------------------------------------//
void DTK_destroyMap( DTK_MapHandle handle )
{
    if ( DataTransferKit::valid_map_handles.count( handle ) )
    {
        auto dtk = reinterpret_cast<DataTransferKit::DTK_Map *>( handle );
        delete dtk;
        DataTransferKit::valid_map_handles.erase( handle );
        errno = DTK_SUCCESS;
    }
    else
    {
        errno = DTK_INVALID_HANDLE;
    }
}

//---------------------------------------------------------------------------//

} // end extern "C"
