//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/**
 * \brief DTK_Mesh_def.hpp
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESH_DEF_HPP
#define DTK_MESH_DEF_HPP

#include "DTK_ConfigDefs.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
template <typename SC, typename LO, typename GO, typename NO>
Mesh<SC, LO, GO, NO>::Mesh( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                            const global_id_view node_ids,
                            const connectivity_view connectivity,
                            const coordinate_view coordinates,
                            const std::vector<shards::CellTopology> &topology )
{
    using connectivity_view_2d =
        Kokkos::View<local_ordinal_type **, device_type>;

    std::set<shards::CellTopology> topo_done;
    const size_t n_cells = topology.size();
    size_t cells_left = n_cells;
    _mesh_blocks = Teuchos::rcp( new std::vector<MeshBlock<SC, LO, GO, NO>>() );
    while ( cells_left > 0 )
    {
        // Search for a topology that we haven't treated yet
        size_t topo_pos = 0;
        while ( topo_done.count( topology[topo_pos] ) == 1 )
            ++topo_pos;
        const shards::CellTopology &current_topo = topology[topo_pos];
        const unsigned int n_nodes_per_cell = current_topo.getNodeCount();

        // Count the number of cells that have the current topology
        size_t n_cells_current_topo = 0;
        for ( size_t i = 0; i < n_cells; ++i )
            if ( topology[i] == current_topo )
                ++n_cells_current_topo;

        // Create the MeshBlock for the current_topology
        size_t pos = 0;
        size_t cell = 0;
        // Transform the 1D view of the connectivity to a 2D view
        connectivity_view_2d connectivity_2d(
            "connectivity_2d", n_cells_current_topo, n_nodes_per_cell );
        for ( size_t i = 0; i < n_cells; ++i )
        {
            // Compare topology using the unique key associated to them.
            if ( topology[i].getKey() == current_topo.getKey() )
            {
                for ( unsigned int n = 0; n < n_nodes_per_cell; ++n )
                    connectivity_2d( cell, n ) = connectivity( pos + n );
                ++cell;
                --cells_left;
                if ( cell == n_cells_current_topo )
                    break;
            }

            pos += topology[i].getNodeCount();
        }

        // Create the MeshBlock
        _mesh_blocks->push_back( MeshBlock<SC, LO, GO, NO>(
            comm, node_ids, connectivity_2d, coordinates, current_topo ) );
        topo_done.insert( current_topo );
    }
}

template <typename SC, typename LO, typename GO, typename NO>
Mesh<SC, LO, GO, NO>::Mesh(
    const Teuchos::RCP<std::vector<MeshBlock<SC, LO, GO, NO>>> mesh_blocks )
    : _mesh_blocks( mesh_blocks )
{
}

template <typename SC, typename LO, typename GO, typename NO>
const Teuchos::RCP<std::vector<MeshBlock<SC, LO, GO, NO>>>
Mesh<SC, LO, GO, NO>::meshBlocks() const
{
    return _mesh_blocks;
}
}

#define DTK_MESH_INSTANT( SCALAR, LO, GO, NODE )                               \
    template class Mesh<SCALAR, LO, GO, NODE>;

#endif
