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
 * \brief DTK supported cell topology types.
 */
#ifndef DTK_CELLTYPES_H
#define DTK_CELLTYPES_H

/*!
 * \brief Cell topology enumeration.
 *
 * The following standard cell topology types are supported by DTK.
 */
typedef enum {
    DTK_TRI_3 = 0,
    DTK_TRI_6,
    DTK_QUAD_4,
    DTK_QUAD_9,
    DTK_TET_4,
    DTK_TET_10,
    DTK_TET_11,
    DTK_HEX_8,
    DTK_HEX_20,
    DTK_HEX_27,
    DTK_PYRAMID_5,
    DTK_PYRAMID_13,
    DTK_WEDGE_6,
    DTK_WEDGE_15,
    DTK_WEDGE_18,
    DTK_N_TOPO
} DTK_CellTopology;

#endif // DTK_CELLTYPES_H
