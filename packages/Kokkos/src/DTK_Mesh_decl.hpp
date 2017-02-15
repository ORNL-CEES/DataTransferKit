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
 * \brief DTK_Mesh_decl.hpp
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESH_DECL_HPP
#define DTK_MESH_DECL_HPP

#include "DTK_ConfigDefs.hpp"
#include "DTK_MeshBlock.hpp"

namespace DataTransferKit
{
template <typename SC, typename LO, typename GO, typename NO>
class Mesh
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
     * Connectivity view. Dimensions: (Point)
     */
    using connectivity_view =
        Kokkos::View<const local_ordinal_type *, device_type>;

    Mesh( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
          const global_id_view node_ids, const connectivity_view connectivity,
          const coordinate_view coordinates,
          const std::vector<shards::CellTopology> &topologies );

    Mesh( Teuchos::RCP<std::vector<MeshBlock<SC, LO, GO, NO>>> mesh_blocks );

    const Teuchos::RCP<std::vector<MeshBlock<SC, LO, GO, NO>>>
    meshBlocks() const;

  private:
    Teuchos::RCP<std::vector<MeshBlock<SC, LO, GO, NO>>> _mesh_blocks;
};
}

#endif
