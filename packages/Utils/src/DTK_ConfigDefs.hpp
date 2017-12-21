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

#include "DataTransferKit_config.hpp"

#ifdef HAVE_DTK_BOOST
#include <boost/current_function.hpp>
#endif

#include <string>

namespace DataTransferKit
{

#include "DTK_Types.h"

// clang-format off
#ifdef HAVE_DTK_BOOST
#define REGION_NAME(x) BOOST_CURRENT_FUNCTION+std::string(":")+std::string(x)
#else
#define REGION_NAME(x) std::string(x)
#endif
// clang-format on

} // end namespace DataTransferKit

#endif // #ifndef DTK_CONFIGDEFS_HPP
