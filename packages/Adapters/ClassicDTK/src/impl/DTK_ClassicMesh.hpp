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
 * \file DTK_ClassicMesh.hpp
 * \author Stuart R. Slattery
 * \brief Mesh manager declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CLASSICMESH_HPP
#define DTK_CLASSICMESH_HPP

#include <unordered_map>

#include "DTK_Types.hpp"
#include "DTK_MeshManager.hpp"
#include "DTK_MeshTraits.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Tuple.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_FieldContainer.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 \class ClassicMesh
 \brief Wrapper for the classic mesh manager.
 */
//---------------------------------------------------------------------------//
template<class Mesh>
class ClassicMesh
{
  public:

    // Typedefs.
    typedef MeshTraits<Mesh> MT;
    typedef typename MT::global_ordinal_type GlobalOrdinal;

    // Constructor.
    ClassicMesh(
        const Teuchos::RCP<MeshManager<Mesh> >& mesh_manager );

    //! Get the number of mesh blocks.
    int getNumBlocks() const
    { return d_mesh_manager->getNumBlocks(); }

    // Get the local number of elements in the mesh.
    GlobalOrdinal localNumElements() const
    { return d_mesh_manager->localNumElements(); }

    // Get the global number of elements in the mesh.
    GlobalOrdinal globalNumElements() const
    { return d_mesh_manager->globalNumElements(); }

    //! Get a block of mesh.
    Teuchos::RCP<Mesh> getBlock( const int block_id ) const
    { return d_mesh_manager->getBlock(block_id); }

    //! Get the communicator for the mesh.
    Teuchos::RCP<const Teuchos::Comm<int> > comm() const
    { return d_mesh_manager->comm(); }

    //! Get the shards topology of a block.
    shards::CellTopology getBlockTopology( const int block_id ) const
    { return d_block_topo[block_id]; }

    //! Get the physical dimension of the mesh.
    int dim() const
    { return d_mesh_manager->dim(); }

    // Given an element id get its block id.
    int elementBlockId( const GlobalOrdinal gid ) const;

    // Given a vertex global id and block id, get the local id in the block.
    int vertexLocalId( const GlobalOrdinal gid, const int block_id ) const;

    // Given a element global id and block id, get the local id in the block.
    int elementLocalId( const GlobalOrdinal gid, const int block_id ) const;

    // Given an element global id and its block id get the coordinates of the
    // element nodes.
    Intrepid::FieldContainer<double>
    getElementNodeCoordinates( const GlobalOrdinal gid,
                               const int block_id ) const;

    // Get the connectivity of an element.
    Teuchos::Array<SupportId>
    getElementConnectivity( const GlobalOrdinal gid, const int block_id ) const;

    // Get the coordinates of a single node.
    Teuchos::Array<double>
    getNodeCoordinates( const GlobalOrdinal gid, const int block_id ) const;

  public:

    // Pointer to element global ids. Public for iterator access.
    Teuchos::Array<Teuchos::ArrayRCP<const GlobalOrdinal> > d_element_gids;

  private:

    // Given a block id create its shards topology.
    shards::CellTopology createBlockTopology( const int block_id ) const;

  private:

    // Classic mesh manager.
    Teuchos::RCP<MeshManager<Mesh> > d_mesh_manager;

    // Block topologies.
    Teuchos::Array<shards::CellTopology> d_block_topo;

    // Element global id to block map.
    std::unordered_map<GlobalOrdinal,int> d_element_block_map;

    // Block-wise vertex global-to-local map.
    Teuchos::Array<std::unordered_map<GlobalOrdinal,int> > d_vertex_g2l;

    // Block-wise element global-to-local map.
    Teuchos::Array<std::unordered_map<GlobalOrdinal,int> > d_element_g2l;

    // Pointer to vertex global ids.
    Teuchos::Array<Teuchos::ArrayRCP<const GlobalOrdinal> > d_vertex_gids;

    // Pointer to vertex coordinates.
    Teuchos::Array<Teuchos::ArrayRCP<const double> > d_vertex_coords;

    // Pointer to element connectivity.
    Teuchos::Array<Teuchos::ArrayRCP<const GlobalOrdinal> > d_element_conn;

    // Pointer to permutation lists.
    Teuchos::Array<Teuchos::ArrayRCP<const int> > d_permutation;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_ClassicMesh_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_CLASSICMESH_HPP

//---------------------------------------------------------------------------//
// end DTK_ClassicMesh.hpp
//---------------------------------------------------------------------------//

