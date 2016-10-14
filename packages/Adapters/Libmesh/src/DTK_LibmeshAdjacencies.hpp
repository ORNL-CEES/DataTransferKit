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
 * \brief DTK_LibmeshAdjacencies.hpp
 * \author Stuart R. Slattery
 * \brief Adjacency data for Libmesh mesh.
 */
//---------------------------------------------------------------------------//

#ifndef LIBMESHDTKADAPTERS_ADJACENCIES_HPP
#define LIBMESHDTKADAPTERS_ADJACENCIES_HPP

#include <unordered_map>

#include <DTK_Types.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_RCP.hpp>

#include <libmesh/elem.h>
#include <libmesh/mesh_base.h>
#include <libmesh/node.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class LibmeshAdjacencies
  \brief Libmesh adjacency information.

  Libmesh doesn't create the upward adjacency graph so this class takes care
  of that. For now, only element and node adjacencies are supported in this
  implementation. This will be sufficient for shared-domain coupling. For
  surface transfers, we will need to add support for edges and faces.
*/
//---------------------------------------------------------------------------//
class LibmeshAdjacencies
{
  public:
    // Constructor.
    LibmeshAdjacencies( const Teuchos::RCP<libMesh::MeshBase> &mesh );

    // Get the adjacency of a libmesh geom object.
    template <class FromGeomType, class ToGeomType>
    void getLibmeshAdjacencies(
        const Teuchos::Ptr<FromGeomType> &entity,
        Teuchos::Array<Teuchos::Ptr<ToGeomType>> &adjacent_entities ) const;

    // Given a node global id get its pointer.
    libMesh::Node *getNodeById( const DataTransferKit::EntityId id ) const;

    // Given a elem global id get its pointer.
    libMesh::Elem *getElemById( const DataTransferKit::EntityId id ) const;

  private:
    // libMesh mesh.
    Teuchos::RCP<libMesh::MeshBase> d_mesh;

    // Node-to-element map.
    std::unordered_multimap<libMesh::Node *, libMesh::Elem *>
        d_node_to_elem_map;

    // Id-to-node map.
    std::unordered_map<DataTransferKit::EntityId, libMesh::Node *>
        d_node_id_map;

    // Id-to-elem map.
    std::unordered_map<DataTransferKit::EntityId, libMesh::Elem *>
        d_elem_id_map;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end LIBMESHDTKADAPTERS_ADJACENCIES_HPP

//---------------------------------------------------------------------------//
// end DTK_LibmeshAdjacencies.hpp
//---------------------------------------------------------------------------//
