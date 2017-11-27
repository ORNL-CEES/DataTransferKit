/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#ifndef DTK_TOPOLOGY_HPP
#define DTK_TOPOLOGY_HPP

#include <DTK_CellTypes.h>

namespace DataTransferKit
{
struct Topology
{
  public:
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
}

#endif
