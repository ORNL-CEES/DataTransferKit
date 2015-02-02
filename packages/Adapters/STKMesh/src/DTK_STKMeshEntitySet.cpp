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
#include "DTK_STKMeshEntity.hpp"
#include "DTK_STKMeshEntityIterator.hpp"
#include "DTK_STKMeshEntityIteratorRange.hpp"
#include "DTK_STKMeshHelpers.hpp"
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
    stk::mesh::Entity stk_entity = 
	d_bulk_data->get_entity( 
	    STKMeshHelpers::getRankFromType(entity_type,physicalDimension()), 
	    entity_id );
    entity = STKMeshEntity( stk_entity,	d_bulk_data.ptr() );
}

//---------------------------------------------------------------------------//
// Get an iterator over a subset of the entity set that satisfies the given
// predicate. 
EntityIterator STKMeshEntitySet::entityIterator(
    const EntityType entity_type,
    const std::function<bool(Entity)>& predicate ) const
{
    stk::mesh::EntityRank rank = 
	STKMeshHelpers::getRankFromType( entity_type, physicalDimension() );
    Teuchos::RCP<STKMeshEntityIteratorRange> iterator_range =
	Teuchos::rcp( new STKMeshEntityIteratorRange() );
    stk::mesh::get_entities( *d_bulk_data, rank, iterator_range->d_stk_entities );
    return STKMeshEntityIterator( iterator_range, d_bulk_data, predicate );
}

//---------------------------------------------------------------------------//
// Given an entity, get the entities of the given type that are adjacent to
// it. 
void STKMeshEntitySet::getAdjacentEntities(
    const Entity& entity,
    const EntityType entity_type,
    Teuchos::Array<Entity>& adjacent_entities ) const
{
    const stk::mesh::Entity& stk_entity = STKMeshHelpers::extractEntity(entity);
    stk::mesh::EntityRank rank = 
	STKMeshHelpers::getRankFromType( entity_type, physicalDimension() );
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

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntitySet.cpp
//---------------------------------------------------------------------------//
