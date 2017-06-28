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
 * \brief SWIG macros.
 */
//---------------------------------------------------------------------------//
#ifndef DTK_SWIG_HPP
#define DTK_SWIG_HPP

#ifndef SWIGEXPORT
#if defined( __GNUC__ )
#define SWIGEXPORT __attribute__( ( visibility( "default" ) ) )
#else
#define SWIGEXPORT
#endif
#endif

/* Intel's compiler complains if a variable which was never initialised is
 * cast to void, which is a common idiom which we use to indicate that we
 * are aware a variable isn't used.  So we just silence that warning.
 * See: https://github.com/swig/swig/issues/192 for more discussion.
 */
#ifdef __INTEL_COMPILER
#pragma warning disable 592
#endif

#endif // end DTK_SWIG_HPP
