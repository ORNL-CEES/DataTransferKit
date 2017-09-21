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
 * \brief DTK interface enumerations.
 */
#ifndef DTK_ENUMERATIONS_H
#define DTK_ENUMERATIONS_H

/*---------------------------------------------------------------------------*/
/*!
 * \brief Cell topology enumeration.
 *
 * The following standard cell topology types are supported by DTK.
 */
typedef enum
{
    DTK_TRI_3,
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
    DTK_WEDGE_18
} DTK_CellTopology;

/*---------------------------------------------------------------------------*/
/*!
 * \brief Discretization type enumeration.
 *
 * The following standard discretization types are supported by DTK.
 */
typedef enum
{
    DTK_NODE,
    DTK_EDGE,
    DTK_FACE,
    DTK_CELL,
    DTK_HGRAD,
    DTK_HDIV,
    DTK_HCURL
} DTK_Discretization;

/*---------------------------------------------------------------------------*/
/*!
 * \brief Cell basis point type enumeration.
 *
 * The following basis point types are supported by DTK for higher order
 * elements.
 */
typedef enum
{
    DTK_EQUISPACED,
    DTK_SPECTRAL,
    DTK_SPECTRAL_OPEN,
    DTK_WARPBLEND
} DTK_BasisPointType;

/*---------------------------------------------------------------------------*/

#endif // DTK_ENUMERATIONS_H
