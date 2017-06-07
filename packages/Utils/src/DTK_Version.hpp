/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#ifndef DTK_VERSION_HPP
#define DTK_VERSION_HPP

#include <DataTransferKit_config.hpp>

#include <string>

namespace DataTransferKit
{

inline std::string version() { return DataTransferKit_VERSION_STRING; }
}

#endif
