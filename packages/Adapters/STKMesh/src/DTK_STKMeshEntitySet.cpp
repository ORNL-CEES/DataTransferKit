//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \brief DTK_STKMeshEntitySet.cpp
 * \author Stuart R. Slattery
 * \brief STK mesh entity set.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "DTK_STKMeshEntitySet.hpp"
#include "DTK_STKMeshEntityExtraData.hpp"
#include "DTK_STKMeshEntityIterator.hpp"
#include "DTK_DBC.hpp"

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <Teuchos_DefaultMpiComm.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
STKMeshEntitySet::STKMeshEntitySet( 
    const Teuchos::RCP<stk::mesh::BulkData>& bulk_data )
    : d_bulk_data( bulk_data )
{ 
    DTK_REQUIRE( Teuchos::nonnull(d_bulk_data) );
}

//---------------------------------------------------------------------------//
// Destructor.
STKMeshEntitySet::~STKMeshEntitySet()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the parallel communicator for the entity set.
Teuchos::RCP<const Teuchos::Comm<int> > STKMeshEntitySet::communicator() const
{
    return Teuchos::rcp( new Teuchos::MpiComm<int>(d_bulk_data->parallel()) );
}

//---------------------------------------------------------------------------//
// Return the largest physical dimension of the entities in the set. 
int STKMeshEntitySet::physicalDimension() const
{
    return d_bulk_data->mesh_meta_data().spatial_dimension();
}

//---------------------------------------------------------------------------//
// Given an EntityId, get the entity.
void STKMeshEntitySet::getEntity( const EntityType entity_type,
				  const EntityId entity_id, 
				  Entity& entity ) const
{
    entity = STKMeshEntity( 
	d_bulk_data->get_entity( getRankFromEntityType(entity_type), entity_id ),
	d_bulk_data.ptr()
	);
}

//---------------------------------------------------------------------------//
// Get an iterator over a subset of the entity set that satisfies the given
// predicate. 
EntityIterator STKMeshEntitySet::entityIterator(
    const EntityType entity_type,
    const std::function<bool(Entity)>& predicate ) const
{
    stk::mesh::EntityRank rank = getRankFromEntityType( entity_type );
    std::vector<stk::mesh::Entity> stk_entities;
    stk::mesh::get_entities( *d_bulk_data, rank, stk_entities );
    return STKMeshEntityIterator( stk_entities, d_bulk_data, predicate );
}

//---------------------------------------------------------------------------//
// Given an entity, get the entities of the given type that are adjacent to
// it. 
void STKMeshEntitySet::getAdjacentEntities(
    const Entity& entity,
    const EntityType entity_type,
    Teuchos::Array<Entity>& adjacent_entities ) const
{
    Teuchos::RCP<EntityExtraData> extra_data = entity.extraData();
    stk::mesh::Entity stk_entity = 
	Teuchos::rcp_dynamic_cast<STKMeshEntityExtraData>(
	    extra_data)->d_stk_entity;
    stk::mesh::EntityRank rank = getRankFromEntityType( entity_type );
    const stk::mesh::Entity* begin = 
	d_bulk_data->begin( stk_entity, rank );
    const stk::mesh::Entity* end = d_bulk_data->end( stk_entity, rank );
    Teuchos::Array<stk::mesh::Entity> stk_adjacencies( begin, end );
    adjacent_entities.resize( stk_adjacencies.size() );
    Teuchos::Array<Entity>::iterator entity_it;
    Teuchos::Array<stk::mesh::Entity>::iterator stk_it;
    for ( entity_it = adjacent_entities.begin(),
	     stk_it = stk_adjacencies.begin();
	  entity_it != adjacent_entities.end();
	  ++entity_it, ++stk_it )
    {
	*entity_it = STKMeshEntity( *stk_it, d_bulk_data.ptr() );
    }
}

//---------------------------------------------------------------------------//
// Given an entity type, get the STK entity rank.
stk::mesh::EntityRank 
STKMeshEntitySet::getRankFromEntityType( const EntityType entity_type ) const
{
    stk::mesh::EntityRank rank = stk::topology::INVALID_RANK;

    switch( entity_type )
    {
	case ENTITY_TYPE_NODE:
	    rank = stk::topology::NODE_RANK;
	    break;
	case ENTITY_TYPE_EDGE:
	    rank = stk::topology::EDGE_RANK;
	    break;
	case ENTITY_TYPE_FACE:
	    rank = stk::topology::FACE_RANK;
	    break;
	case ENTITY_TYPE_VOLUME:
	    rank = stk::topology::ELEM_RANK;
	    break;
	default:
	    DTK_CHECK( ENTITY_TYPE_NODE == entity_type ||
		       ENTITY_TYPE_EDGE == entity_type ||
		       ENTITY_TYPE_FACE == entity_type ||
		       ENTITY_TYPE_VOLUME == entity_type );
	    break;
    }

    return rank;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntitySet.cpp
//---------------------------------------------------------------------------//
