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

#ifndef DTK_TOPOLOGY_HPP
#define DTK_TOPOLOGY_HPP

#include <DTK_CellTypes.h>

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

#include <array>

namespace DataTransferKit
{
struct Topology
{
    Topology() = default;

    Topology( unsigned int d, unsigned int n, DTK_CellTopology t )
        : dim( d )
        , n_nodes( n )
        , topo( t )
    {
    }

    unsigned int dim;
    unsigned int n_nodes;
    DTK_CellTopology topo;
};

class Topologies
{
  public:
    Topologies()
    {
        _topologies[DTK_TRI_3] = Topology( 2, 3, DTK_TRI_3 );
        _topologies[DTK_TRI_6] = Topology( 2, 6, DTK_TRI_6 );
        _topologies[DTK_QUAD_4] = Topology( 2, 4, DTK_QUAD_4 );
        _topologies[DTK_QUAD_9] = Topology( 2, 9, DTK_QUAD_9 );
        _topologies[DTK_TET_4] = Topology( 3, 4, DTK_TET_4 );
        _topologies[DTK_TET_10] = Topology( 3, 10, DTK_TET_10 );
        _topologies[DTK_TET_11] = Topology( 3, 11, DTK_TET_11 );
        _topologies[DTK_HEX_8] = Topology( 3, 8, DTK_HEX_8 );
        _topologies[DTK_HEX_20] = Topology( 3, 20, DTK_HEX_20 );
        _topologies[DTK_HEX_27] = Topology( 3, 27, DTK_HEX_27 );
        _topologies[DTK_PYRAMID_5] = Topology( 3, 5, DTK_PYRAMID_5 );
        _topologies[DTK_PYRAMID_13] = Topology( 3, 13, DTK_PYRAMID_13 );
        _topologies[DTK_WEDGE_6] = Topology( 3, 6, DTK_WEDGE_6 );
        _topologies[DTK_WEDGE_15] = Topology( 3, 15, DTK_WEDGE_15 );
        _topologies[DTK_WEDGE_18] = Topology( 3, 18, DTK_WEDGE_18 );
    }

    Topology &operator[]( int const i ) { return _topologies[i]; }
    Topology const &operator[]( int const i ) const { return _topologies[i]; }

  private:
    std::array<Topology, DTK_N_TOPO> _topologies;
};

struct HEX_8
{
    typedef Intrepid2::Impl::Basis_HGRAD_HEX_C1_FEM basis_type;
    typedef Intrepid2::Impl::Hexahedron<8> topo_type;
};

struct HEX_27
{
    typedef Intrepid2::Impl::Basis_HGRAD_HEX_DEG2_FEM<false> basis_type;
    typedef Intrepid2::Impl::Hexahedron<27> topo_type;
};

struct PYRAMID_5
{
    typedef Intrepid2::Impl::Basis_HGRAD_PYR_C1_FEM basis_type;
    typedef Intrepid2::Impl::Pyramid<5> topo_type;
};

struct QUAD_4
{
    typedef Intrepid2::Impl::Basis_HGRAD_QUAD_C1_FEM basis_type;
    typedef Intrepid2::Impl::Quadrilateral<4> topo_type;
};

struct QUAD_9
{
    typedef Intrepid2::Impl::Basis_HGRAD_QUAD_DEG2_FEM<false> basis_type;
    typedef Intrepid2::Impl::Quadrilateral<9> topo_type;
};

struct TET_4
{
    typedef Intrepid2::Impl::Basis_HGRAD_TET_C1_FEM basis_type;
    typedef Intrepid2::Impl::Tetrahedron<4> topo_type;
};

struct TET_10
{
    typedef Intrepid2::Impl::Basis_HGRAD_TET_C2_FEM basis_type;
    typedef Intrepid2::Impl::Tetrahedron<10> topo_type;
};

struct TRI_3
{
    typedef Intrepid2::Impl::Basis_HGRAD_TRI_C1_FEM basis_type;
    typedef Intrepid2::Impl::Triangle<3> topo_type;
};

struct TRI_6
{
    typedef Intrepid2::Impl::Basis_HGRAD_TRI_C2_FEM basis_type;
    typedef Intrepid2::Impl::Triangle<6> topo_type;
};

struct WEDGE_6
{
    typedef Intrepid2::Impl::Basis_HGRAD_WEDGE_C1_FEM basis_type;
    typedef Intrepid2::Impl::Wedge<6> topo_type;
};

struct WEDGE_18
{
    typedef Intrepid2::Impl::Basis_HGRAD_WEDGE_DEG2_FEM<false> basis_type;
    typedef Intrepid2::Impl::Wedge<18> topo_type;
};
} // namespace DataTransferKit

#endif
