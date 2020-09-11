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

#ifndef DTK_VERSION_HPP
#define DTK_VERSION_HPP

#include <DTK_ConfigDefs.hpp>

#include <string>

namespace DataTransferKit
{

inline std::string version() { return DataTransferKit_VERSION_STRING; }

inline std::string gitCommitHash() { return DataTransferKit_GIT_COMMIT_HASH; }
} // namespace DataTransferKit

#endif
