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

#include <DTK_FE.hpp>

namespace DataTransferKit
{
FE getFE( DTK_CellTopology topo, DTK_FEType fe_type )
{
    if ( topo == DTK_HEX_8 )
    {
        switch ( fe_type )
        {
        case DTK_HCURL:
        {
            return FE::HEX_HCURL_1;
        }
        case DTK_HDIV:
        {
            return FE::HEX_HDIV_1;
        }
        case DTK_HGRAD:
        {
            return FE::HEX_HGRAD_1;
        }
        default:
            return FE::DUMMY;
        }
    }
    else if ( topo == DTK_HEX_27 )
    {
        if ( fe_type == DTK_HGRAD )
            return FE::HEX_HGRAD_2;
        else
            return FE::DUMMY;
    }
    else if ( ( topo == DTK_PYRAMID_5 ) && ( fe_type == DTK_HGRAD ) )
    {
        return FE::PYR_HGRAD_1;
    }
    else if ( topo == DTK_QUAD_4 )
    {
        switch ( fe_type )
        {
        case DTK_HCURL:
        {
            return FE::QUAD_HCURL_1;
        }
        case DTK_HDIV:
        {
            return FE::QUAD_HDIV_1;
        }
        case DTK_HGRAD:
        {
            return FE::QUAD_HGRAD_1;
        }
        default:
            return FE::DUMMY;
        }
    }
    else if ( ( topo == DTK_QUAD_9 ) && ( fe_type == DTK_HGRAD ) )
    {
        return FE::QUAD_HGRAD_2;
    }
    else if ( topo == DTK_TET_4 )
    {
        switch ( fe_type )
        {
        case DTK_HCURL:
        {
            return FE::TET_HCURL_1;
        }
        case DTK_HDIV:
        {
            return FE::TET_HDIV_1;
        }
        case DTK_HGRAD:
        {
            return FE::TET_HGRAD_1;
        }
        default:
            return FE::DUMMY;
        }
    }
    else if ( ( topo == DTK_TET_10 ) && ( fe_type == DTK_HGRAD ) )
    {
        return FE::TET_HGRAD_2;
    }
    else if ( ( topo == DTK_TRI_3 ) && ( fe_type == DTK_HGRAD ) )
    {
        return FE::TRI_HGRAD_1;
    }
    else if ( ( topo == DTK_TRI_6 ) && ( fe_type == DTK_HGRAD ) )
    {
        return FE::TRI_HGRAD_2;
    }
    else if ( ( topo == DTK_WEDGE_6 ) && ( fe_type == DTK_HGRAD ) )
    {
        return FE::WEDGE_HGRAD_1;
    }
    else if ( ( topo == DTK_WEDGE_18 ) && ( fe_type == DTK_HGRAD ) )
    {
        return FE::WEDGE_HGRAD_2;
    }
    else
        return FE::DUMMY;
}
} // namespace DataTransferKit
