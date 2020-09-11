/****************************************************************************
 * Copyright (c) 2012-2020 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
/*!
 * \file
 * \brief C adapter to UserFunctionRegistry.
 */
#ifndef DTK_C_API_HPP
#define DTK_C_API_HPP

#include <memory>

#include <DataTransferKit_config.hpp>

#include "DTK_C_API.h"
#include "DTK_UserFunctionRegistry.hpp"

namespace DataTransferKit
{

struct DTK_Registry
{
    DTK_Registry( DTK_MemorySpace space )
    {
        _registry = std::make_shared<UserFunctionRegistry<double>>();
        _space = space;
    }

    std::shared_ptr<UserFunctionRegistry<double>> _registry;
    DTK_MemorySpace _space;
};
} // namespace DataTransferKit

#endif // DTK_C_API_HPP
