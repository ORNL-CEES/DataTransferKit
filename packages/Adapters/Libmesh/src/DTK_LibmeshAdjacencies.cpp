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
 * \brief DTK_LibmeshAdjacencies.cpp
 * \author Stuart R. Slattery
 * \brief Adjacency data for Libmesh mesh.
 */
//---------------------------------------------------------------------------//

#include "DTK_LibmeshAdjacencies.hpp"
#include <DTK_DBC.hpp>

#include <libmesh/edge.h>
#include <libmesh/face.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
LibmeshAdjacencies::LibmeshAdjacencies(
    const Teuchos::RCP<libMesh::MeshBase> &mesh )
    : d_mesh( mesh )
{
    // Map nodes to elements and elements to their ids.
    int num_nodes = 0;
    libMesh::MeshBase::element_iterator elem_begin =
        d_mesh->local_elements_begin();
    libMesh::MeshBase::element_iterator elem_end = d_mesh->local_elements_end();
    for ( auto elem = elem_begin; elem != elem_end; ++elem )
    {
        d_elem_id_map.emplace( ( *elem )->id(), *elem );
        num_nodes = ( *elem )->n_nodes();
        for ( int n = 0; n < num_nodes; ++n )
        {
            d_node_to_elem_map.emplace( ( *elem )->get_node( n ), *elem );
        }
    }

    // Map nodes to their ids.
    libMesh::MeshBase::node_iterator node_begin = d_mesh->local_nodes_begin();
    libMesh::MeshBase::node_iterator node_end = d_mesh->local_nodes_end();
    for ( auto node = node_begin; node != node_end; ++node )
    {
        d_node_id_map.emplace( ( *node )->id(), *node );
    }
}

//---------------------------------------------------------------------------//
/*!
 * Get the adjacency of a libmesh geom object. Node to elem overload.
 */
template <>
void LibmeshAdjacencies::getLibmeshAdjacencies<libMesh::Node, libMesh::Elem>(
    const Teuchos::Ptr<libMesh::Node> &entity,
    Teuchos::Array<Teuchos::Ptr<libMesh::Elem>> &adjacent_entities ) const
{
    auto elem_range = d_node_to_elem_map.equal_range( entity.getRawPtr() );
    int num_elem = std::distance( elem_range.first, elem_range.second );
    adjacent_entities.resize( num_elem );
    int e = 0;
    for ( auto node_elems = elem_range.first; node_elems != elem_range.second;
          ++node_elems, ++e )
    {
        adjacent_entities[e] = Teuchos::ptr( node_elems->second );
    }
}

//---------------------------------------------------------------------------//
/*!
 * Get the adjacency of a libmesh geom object. Elem to node overload.
 */
template <>
void LibmeshAdjacencies::getLibmeshAdjacencies<libMesh::Elem, libMesh::Node>(
    const Teuchos::Ptr<libMesh::Elem> &entity,
    Teuchos::Array<Teuchos::Ptr<libMesh::Node>> &adjacent_entities ) const
{
    int num_nodes = entity->n_nodes();
    adjacent_entities.resize( num_nodes );
    for ( int n = 0; n < num_nodes; ++n )
    {
        adjacent_entities[n] = Teuchos::ptr( entity->get_node( n ) );
    }
}

//---------------------------------------------------------------------------//
/*!
 * Get the adjacency of a libmesh geom object. Node to node overload.
 */
template <>
void LibmeshAdjacencies::getLibmeshAdjacencies<libMesh::Node, libMesh::Node>(
    const Teuchos::Ptr<libMesh::Node> & /*entity*/,
    Teuchos::Array<Teuchos::Ptr<libMesh::Node>> &adjacent_entities ) const
{
    adjacent_entities.clear();
}

//---------------------------------------------------------------------------//
/*!
 * Get the adjacency of a libmesh geom object. Elem to elem overload.
 */
template <>
void LibmeshAdjacencies::getLibmeshAdjacencies<libMesh::Elem, libMesh::Elem>(
    const Teuchos::Ptr<libMesh::Elem> & /*entity*/,
    Teuchos::Array<Teuchos::Ptr<libMesh::Elem>> &adjacent_entities ) const
{
    adjacent_entities.clear();
}

//---------------------------------------------------------------------------//
// Given a node global id get its pointer.
libMesh::Node *
LibmeshAdjacencies::getNodeById( const DataTransferKit::EntityId id ) const
{
    DTK_REQUIRE( d_node_id_map.count( id ) );
    return d_node_id_map.find( id )->second;
}

//---------------------------------------------------------------------------//
// Given a elem global id get its pointer.
libMesh::Elem *
LibmeshAdjacencies::getElemById( const DataTransferKit::EntityId id ) const
{
    DTK_REQUIRE( d_elem_id_map.count( id ) );
    return d_elem_id_map.find( id )->second;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_LibmeshAdjacencies.cpp
//---------------------------------------------------------------------------//
