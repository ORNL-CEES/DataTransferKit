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
/*!
 * \brief DTK_MeshBlock_def.hpp
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHBLOCK_DEF_HPP
#define DTK_MESHBLOCK_DEF_HPP

#include "DTK_ConfigDefs.hpp"

namespace DataTransferKit
{
template <typename SC, typename LO, typename GO, typename NO>
MeshBlock<SC, LO, GO, NO>::MeshBlock(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
    const global_id_view node_ids, const connectivity_view connectivity,
    const coordinate_view coordinates, const shards::CellTopology &topology )
    : _comm( comm )
    , _node_ids( node_ids )
    , _connectivity( connectivity )
    , _coordinates( coordinates )
    , _topology( topology )
{
}

template <typename SC, typename LO, typename GO, typename NO>
size_t MeshBlock<SC, LO, GO, NO>::spaceDim() const
{
    return _coordinates.extent( 1 );
}

template <typename SC, typename LO, typename GO, typename NO>
size_t MeshBlock<SC, LO, GO, NO>::numLocalCells() const
{
    return _connectivity.extent( 0 );
}

template <typename SC, typename LO, typename GO, typename NO>
global_size_t MeshBlock<SC, LO, GO, NO>::numGlobalCells() const
{
    size_t local_num = numLocalCells();
    global_size_t global_num = 0;
    Teuchos::reduceAll( *_comm, Teuchos::REDUCE_SUM, local_num,
                        Teuchos::ptrFromRef( global_num ) );
    return global_num;
}

template <typename SC, typename LO, typename GO, typename NO>
size_t MeshBlock<SC, LO, GO, NO>::numLocalNodes() const
{
    return _coordinates.extent( 0 );
}

template <typename SC, typename LO, typename GO, typename NO>
global_size_t MeshBlock<SC, LO, GO, NO>::numGlobalNodes() const
{
    size_t local_num = numLocalNodes();
    global_size_t global_num = 0;
    Teuchos::reduceAll( *_comm, Teuchos::REDUCE_SUM, local_num,
                        Teuchos::ptrFromRef( global_num ) );
    return global_num;
}

template <typename SC, typename LO, typename GO, typename NO>
const typename MeshBlock<SC, LO, GO, NO>::global_id_view
MeshBlock<SC, LO, GO, NO>::nodeIds() const
{
    return _node_ids;
}

template <typename SC, typename LO, typename GO, typename NO>
const typename MeshBlock<SC, LO, GO, NO>::connectivity_view
MeshBlock<SC, LO, GO, NO>::connectivity() const
{
    return _connectivity;
}

template <typename SC, typename LO, typename GO, typename NO>
const typename MeshBlock<SC, LO, GO, NO>::coordinate_view
MeshBlock<SC, LO, GO, NO>::coordinates() const
{
    return _coordinates;
}

template <typename SC, typename LO, typename GO, typename NO>
const shards::CellTopology &MeshBlock<SC, LO, GO, NO>::topology() const
{
    return _topology;
}
}

#define DTK_MESHBLOCK_INSTANT( SCALAR, LO, GO, NODE )                          \
    template class MeshBlock<SCALAR, LO, GO, NODE>;

#endif
