/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
/*!
 * \file DTK_ConfigDefs.hpp
 * \brief Kokkos helpers.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CONFIGDEFS_HPP
#define DTK_CONFIGDEFS_HPP

#include "DataTransferKitUtils_config.hpp"

#include <boost/current_function.hpp>
#include <cstdint>
#include <string>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//

//! Coordinate typedef.
using Coordinate = double;

//! Local ordinal typedef.
using LocalOrdinal = unsigned int;

//! Global ordinal typedef.
using GlobalOrdinal = uint64_t;

// clang-format off
#define REGION_NAME(x) BOOST_CURRENT_FUNCTION+std::string(":")+std::string(x)
// clang-format on

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_CONFIGDEFS_HPP

//---------------------------------------------------------------------------//
// end DTK_ConfigDefs.hpp
//---------------------------------------------------------------------------//
