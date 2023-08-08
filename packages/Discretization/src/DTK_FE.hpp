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

#ifndef DTK_FE_HPP
#define DTK_FE_HPP

#include <DTK_DBC.hpp>
#include <DTK_FETypes.h>
#include <DTK_Topology.hpp>

#include <Intrepid2_HCURL_HEX_I1_FEM.hpp>
#include <Intrepid2_HCURL_QUAD_I1_FEM.hpp>
#include <Intrepid2_HCURL_TET_I1_FEM.hpp>
#include <Intrepid2_HDIV_HEX_I1_FEM.hpp>
#include <Intrepid2_HDIV_QUAD_I1_FEM.hpp>
#include <Intrepid2_HDIV_TET_I1_FEM.hpp>
#include <Intrepid2_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid2_HGRAD_HEX_C2_FEM.hpp>
#include <Intrepid2_HGRAD_PYR_C1_FEM.hpp>
#include <Intrepid2_HGRAD_QUAD_C1_FEM.hpp>
#include <Intrepid2_HGRAD_QUAD_C2_FEM.hpp>
#include <Intrepid2_HGRAD_TET_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TET_C2_FEM.hpp>
#include <Intrepid2_HGRAD_TRI_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TRI_C2_FEM.hpp>
#include <Intrepid2_HGRAD_WEDGE_C1_FEM.hpp>
#include <Intrepid2_HGRAD_WEDGE_C2_FEM.hpp>

#include <memory>

namespace DataTransferKit
{
enum class FE {
    HEX_HCURL_1,
    HEX_HDIV_1,
    HEX_HGRAD_1,
    HEX_HGRAD_2,
    PYR_HGRAD_1,
    QUAD_HCURL_1,
    QUAD_HDIV_1,
    QUAD_HGRAD_1,
    QUAD_HGRAD_2,
    TET_HCURL_1,
    TET_HDIV_1,
    TET_HGRAD_1,
    TET_HGRAD_2,
    TRI_HGRAD_1,
    TRI_HGRAD_2,
    WEDGE_HGRAD_1,
    WEDGE_HGRAD_2,
    DUMMY
};

struct HEX_HCURL_1
{
    typedef Intrepid2::Impl::Basis_HCURL_HEX_I1_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HCURL_HEX_I1_FEM<T1, T2, T3>;
};

struct HEX_HDIV_1
{
    typedef Intrepid2::Impl::Basis_HDIV_HEX_I1_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HDIV_HEX_I1_FEM<T1, T2, T3>;
};

struct HEX_HGRAD_1
{
    typedef Intrepid2::Impl::Basis_HGRAD_HEX_C1_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HGRAD_HEX_C1_FEM<T1, T2, T3>;
};

struct HEX_HGRAD_2
{
    typedef Intrepid2::Impl::Basis_HGRAD_HEX_DEG2_FEM<false>::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HGRAD_HEX_C2_FEM<T1, T2, T3>;
};

struct PYR_HGRAD_1
{
    typedef Intrepid2::Impl::Basis_HGRAD_PYR_C1_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HGRAD_PYR_C1_FEM<T1, T2, T3>;
};

struct QUAD_HCURL_1
{
    typedef Intrepid2::Impl::Basis_HCURL_QUAD_I1_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HCURL_QUAD_I1_FEM<T1, T2, T3>;
};

struct QUAD_HDIV_1
{
    typedef Intrepid2::Impl::Basis_HDIV_QUAD_I1_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HDIV_QUAD_I1_FEM<T1, T2, T3>;
};

struct QUAD_HGRAD_1
{
    typedef Intrepid2::Impl::Basis_HGRAD_QUAD_C1_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HGRAD_QUAD_C1_FEM<T1, T2, T3>;
};

struct QUAD_HGRAD_2
{
    typedef Intrepid2::Impl::Basis_HGRAD_QUAD_DEG2_FEM<false>::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HGRAD_QUAD_C2_FEM<T1, T2, T3>;
};

struct TET_HCURL_1
{
    typedef Intrepid2::Impl::Basis_HCURL_TET_I1_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HCURL_TET_I1_FEM<T1, T2, T3>;
};

struct TET_HDIV_1
{
    typedef Intrepid2::Impl::Basis_HDIV_TET_I1_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HDIV_TET_I1_FEM<T1, T2, T3>;
};

struct TET_HGRAD_1
{
    typedef Intrepid2::Impl::Basis_HGRAD_TET_C1_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HGRAD_TET_C1_FEM<T1, T2, T3>;
};

struct TET_HGRAD_2
{
    typedef Intrepid2::Impl::Basis_HGRAD_TET_C2_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HGRAD_TET_C2_FEM<T1, T2, T3>;
};

struct TRI_HGRAD_1
{
    typedef Intrepid2::Impl::Basis_HGRAD_TRI_C1_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HGRAD_TRI_C1_FEM<T1, T2, T3>;
};

struct TRI_HGRAD_2
{
    typedef Intrepid2::Impl::Basis_HGRAD_TRI_C2_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HGRAD_TRI_C2_FEM<T1, T2, T3>;
};

struct WEDGE_HGRAD_1
{
    typedef Intrepid2::Impl::Basis_HGRAD_WEDGE_C1_FEM::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HGRAD_WEDGE_C1_FEM<T1, T2, T3>;
};

struct WEDGE_HGRAD_2
{
    typedef Intrepid2::Impl::Basis_HGRAD_WEDGE_DEG2_FEM<false>::Serial<
        Intrepid2::OPERATOR_VALUE>
        feop_type;

    template <typename T1, typename T2, typename T3>
    using basis_type = Intrepid2::Basis_HGRAD_WEDGE_C2_FEM<T1, T2, T3>;
};

struct DUMMY
{
    typedef void feop_type;
};

FE getFE( DTK_CellTopology topo, DTK_FEType fe_type );

/**
 * Return the number of degrees of freedom per cell for a given
 */
template <typename DeviceType>
unsigned int getCardinality( FE fe )
{
    typedef typename Kokkos::View<double *, DeviceType>::traits::execution_space
        ExecutionSpace;
    std::unique_ptr<Intrepid2::Basis<ExecutionSpace, double, double>> basis;
    switch ( fe )
    {
    case FE::HEX_HCURL_1:
    {
        basis.reset(
            new HEX_HCURL_1::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::HEX_HDIV_1:
    {
        basis.reset(
            new HEX_HDIV_1::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::HEX_HGRAD_1:
    {
        basis.reset(
            new HEX_HGRAD_1::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::HEX_HGRAD_2:
    {
        basis.reset(
            new HEX_HGRAD_2::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::PYR_HGRAD_1:
    {
        basis.reset(
            new PYR_HGRAD_1::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::QUAD_HCURL_1:
    {
        basis.reset(
            new QUAD_HCURL_1::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::QUAD_HDIV_1:
    {
        basis.reset(
            new QUAD_HDIV_1::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::QUAD_HGRAD_1:
    {
        basis.reset(
            new QUAD_HGRAD_1::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::QUAD_HGRAD_2:
    {
        basis.reset(
            new QUAD_HGRAD_2::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::TET_HCURL_1:
    {
        basis.reset(
            new QUAD_HCURL_1::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::TET_HDIV_1:
    {
        basis.reset(
            new TET_HDIV_1::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::TET_HGRAD_1:
    {
        basis.reset(
            new TET_HGRAD_1::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::TET_HGRAD_2:
    {
        basis.reset(
            new TET_HGRAD_2::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::TRI_HGRAD_1:
    {
        basis.reset(
            new TRI_HGRAD_1::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::TRI_HGRAD_2:
    {
        basis.reset(
            new TRI_HGRAD_2::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::WEDGE_HGRAD_1:
    {
        basis.reset(
            new WEDGE_HGRAD_1::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    case FE::WEDGE_HGRAD_2:
    {
        basis.reset(
            new WEDGE_HGRAD_2::basis_type<ExecutionSpace, double, double>() );

        break;
    }
    default:
        return 0;
    }

    return basis->getCardinality();
}
} // namespace DataTransferKit

#endif
