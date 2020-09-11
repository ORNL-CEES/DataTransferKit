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
 * \file  DTK_SanitizerMacros.hpp
 * \author Bruno Turcksin
 * \brief  Macros to suppress sanitizer checks.
 */
//---------------------------------------------------------------------------//

#if defined( __clang__ )
#define IGNORE_UNDEFINED_SANITIZE                                              \
    __attribute__( ( no_sanitize( "undefined" ) ) )
#else
#define IGNORE_UNDEFINED_SANITIZE
#endif
