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
 * \brief DTK_LibmeshEntitySet.cpp
 * \author Stuart R. Slattery
 * \brief Libmesh mesh entity set.
 */
//---------------------------------------------------------------------------//

#include "DTK_LibmeshEntitySet.hpp"
#include "DTK_LibmeshEntityIterator.hpp"

#include <Teuchos_DefaultMpiComm.hpp>

#include <libmesh/elem.h>
#include <libmesh/node.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
LibmeshEntitySet::LibmeshEntitySet(
    const Teuchos::RCP<libMesh::MeshBase> &libmesh_mesh )
    : d_libmesh_mesh( libmesh_mesh )
    , d_adjacencies( new LibmeshAdjacencies( libmesh_mesh ) )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Get the parallel communicator for the entity set.
Teuchos::RCP<const Teuchos::Comm<int>> LibmeshEntitySet::communicator() const
{
    return Teuchos::rcp(
        new Teuchos::MpiComm<int>( d_libmesh_mesh->comm().get() ) );
}

//---------------------------------------------------------------------------//
// Return the largest physical dimension of the entities in the set.
int LibmeshEntitySet::physicalDimension() const
{
    return d_libmesh_mesh->mesh_dimension();
}

//---------------------------------------------------------------------------//
// Given an EntityId, get the entity.
void LibmeshEntitySet::getEntity( const EntityId entity_id,
                                  const int topological_dimension,
                                  Entity &entity ) const
{
    if ( 0 == topological_dimension )
    {
        entity = LibmeshEntity<libMesh::Node>(
            Teuchos::ptr( d_adjacencies->getNodeById( entity_id ) ),
            d_libmesh_mesh.ptr(), d_adjacencies.ptr() );
    }
    else
    {
        entity = LibmeshEntity<libMesh::Elem>(
            Teuchos::ptr( d_adjacencies->getElemById( entity_id ) ),
            d_libmesh_mesh.ptr(), d_adjacencies.ptr() );
    }
}

//---------------------------------------------------------------------------//
// Get an iterator over a subset of the entity set that satisfies the given
// predicate.
EntityIterator LibmeshEntitySet::entityIterator(
    const int topological_dimension,
    const std::function<bool( Entity )> &predicate ) const
{
    EntityIterator entity_it;
    if ( 0 == topological_dimension )
    {
        entity_it =
            LibmeshEntityIterator<libMesh::MeshBase::const_node_iterator>(
                d_libmesh_mesh->local_nodes_begin(),
                d_libmesh_mesh->local_nodes_begin(),
                d_libmesh_mesh->local_nodes_end(), d_libmesh_mesh.ptr(),
                d_adjacencies.ptr(), predicate );
    }
    else
    {
        entity_it =
            LibmeshEntityIterator<libMesh::MeshBase::const_element_iterator>(
                d_libmesh_mesh->local_elements_begin(),
                d_libmesh_mesh->local_elements_begin(),
                d_libmesh_mesh->local_elements_end(), d_libmesh_mesh.ptr(),
                d_adjacencies.ptr(), predicate );
    }
    return entity_it;
}

//---------------------------------------------------------------------------//
// Given an entity, get the entities of the given type that are adjacent to
// it.
void LibmeshEntitySet::getAdjacentEntities(
    const Entity &entity, const int adjacent_dimension,
    Teuchos::Array<Entity> &adjacent_entities ) const
{
    int entity_topo_dim = entity.topologicalDimension();

    if ( 0 == entity_topo_dim )
    {
        if ( 0 == adjacent_dimension )
        {
            adjacent_entities.clear();
        }
        else
        {
            getAdjacentEntitiesImpl<libMesh::Node, libMesh::Elem>(
                entity, adjacent_entities );
        }
    }
    else
    {
        if ( 0 == adjacent_dimension )
        {
            getAdjacentEntitiesImpl<libMesh::Elem, libMesh::Node>(
                entity, adjacent_entities );
        }
        else
        {
            getAdjacentEntitiesImpl<libMesh::Elem, libMesh::Elem>(
                entity, adjacent_entities );
        }
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_LibmeshEntitySet.cpp
//---------------------------------------------------------------------------//
