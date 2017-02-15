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
 * \brief DTK_MeshBlock_decl.hpp
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHBLOCK_DECL_HPP
#define DTK_MESHBLOCK_DECL_HPP

#include "DTK_ConfigDefs.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

#include <Kokkos_Core.hpp>

#include <Shards_CellTopology.hpp>

namespace DataTransferKit
{
template <typename SC, typename LO, typename GO, typename NO>
class MeshBlock
{
  public:
    using scalar_type = SC;
    using local_ordinal_type = LO;
    using global_ordinal_type = GO;
    using node_type = NO;
    using device_type = typename NO::device_type;
    using execution_space = typename device_type::execution_space;
    using memory_space = typename device_type::memory_space;

    /**
     *  Point global ids view. Dimensions: (Point)
     */
    using global_id_view =
        Kokkos::View<const global_ordinal_type *, device_type>;

    /**
     * Point coordinates view. Dimensions: (Point,SpaceDim)
     */
    using coordinate_view = Kokkos::View<const scalar_type **, device_type>;

    /**
     * Connectivity view. Dimensions: (Cell,Point)
     */

    using connectivity_view =
        Kokkos::View<const local_ordinal_type **, device_type>;

    MeshBlock( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
               const global_id_view node_ids,
               const connectivity_view connectivity,
               const coordinate_view coordinates,
               const shards::CellTopology &topology );

    // Get the spatial dimension of the domain.
    size_t spaceDim() const;

    // Get the local number of cells.
    size_t numLocalCells() const;

    // Get the global number of cells.
    global_size_t numGlobalCells() const;

    // Get the local number of nodes.
    size_t numLocalNodes() const;

    // Get the global number of nodes.
    global_size_t numGlobalNodes() const;

    // Get the nodes global ids.
    const global_id_view nodeIds() const;

    const connectivity_view connectivity() const;

    // Get the the point coordinates.
    const coordinate_view coordinates() const;

    const shards::CellTopology &topology() const;

  private:
    Teuchos::RCP<const Teuchos::Comm<int>> _comm;

    const global_id_view _node_ids;

    const connectivity_view _connectivity;

    const coordinate_view _coordinates;

    /**
     * This is necessary because if a cell has 4 nodes, we don't
     * know if it is a quadrilateral or a triangle with a curved edge.
     */
    shards::CellTopology _topology;
};
}

#endif
