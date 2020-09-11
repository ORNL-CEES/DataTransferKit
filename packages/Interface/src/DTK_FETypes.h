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
 * \file
 * \brief DTK supported finite element types.
 */
#ifndef DTK_FETYPES_H
#define DTK_FETYPES_H

/*!
 * \brief Finite element types enumeration.
 *
 * The following standard finite element types are supported by DTK.
 */
typedef enum { DTK_HGRAD = 0, DTK_HDIV, DTK_HCURL, DTK_N_FEM } DTK_FEType;

#endif // DTK_FETYPES_H
