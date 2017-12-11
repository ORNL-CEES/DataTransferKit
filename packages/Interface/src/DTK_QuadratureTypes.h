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
 * \brief DTK supported quadrature types.
 */
#ifndef DTK_QUADRATURETYPES_H
#define DTK_QUADRATURETYPES_H

/*!
 * \brief Quadrature types enumeration.
 *
 * The following standard quadrature types are supported by DTK.
 */
typedef enum {
    DTK_EQUISPACED = 0,
    DTK_WARPBLEND,
    DTK_GAUSS,
    DTK_N_QUAD
} DTK_Quadrature;

#endif // DTK_FETYPES_H
