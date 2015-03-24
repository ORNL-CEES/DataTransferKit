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
 * \brief DTK_MoabEntitySet.cpp
 * \author Stuart R. Slattery
 * \brief Moab mesh entity set.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "DTK_MoabEntitySet.hpp"
#include "DTK_MoabEntity.hpp"
#include "DTK_MoabEntityIterator.hpp"
#include "DTK_MoabEntityIteratorRange.hpp"
#include "DTK_MoabHelpers.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_DefaultMpiComm.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
MoabEntitySet::MoabEntitySet( 
    const Teuchos::RCP<moab::ParallelComm>& moab_mesh )
    : d_moab_mesh( moab_mesh )
    , d_set_indexer( Teuchos::rcp(new MoabMeshSetIndexer(moab_mesh)) )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the parallel communicator for the entity set.
Teuchos::RCP<const Teuchos::Comm<int> > MoabEntitySet::communicator() const
{
    return Teuchos::rcp( new Teuchos::MpiComm<int>(d_moab_mesh->comm()) );
}

//---------------------------------------------------------------------------//
// Return the largest physical dimension of the entities in the set. 
int MoabEntitySet::physicalDimension() const
{
    int dimension = 0;
    DTK_CHECK_ERROR_CODE(
	d_moab_mesh->get_moab()->get_dimension( dimension )
	);
    return dimension;
}

//---------------------------------------------------------------------------//
// Given an EntityId, get the entity.
void MoabEntitySet::getEntity( const EntityType entity_type,
			       const EntityId entity_id, 
			       Entity& entity ) const
{
    entity = MoabEntity( Teuchos::as<moab::EntityHandle>(entity_id),
			 d_moab_mesh.ptr(),
			 d_set_indexer.ptr() );
}

//---------------------------------------------------------------------------//
// Get an iterator over a subset of the entity set that satisfies the given
// predicate. 
EntityIterator MoabEntitySet::entityIterator(
    const EntityType entity_type,
    const PredicateFunction& predicate ) const
{
    Teuchos::RCP<MoabEntityIteratorRange> iterator_range =
	Teuchos::rcp( new MoabEntityIteratorRange() );
    DTK_CHECK_ERROR_CODE(
	d_moab_mesh->get_moab()->get_entities_by_dimension(
	    0,
	    Teuchos::as<int>( entity_type ),
	    iterator_range->d_moab_entities )
	);
    return MoabEntityIterator( 
	iterator_range, d_moab_mesh.ptr(), d_set_indexer.ptr(), predicate );
}

//---------------------------------------------------------------------------//
// Given an entity, get the entities of the given type that are adjacent to
// it. 
void MoabEntitySet::getAdjacentEntities(
    const Entity& entity,
    const EntityType entity_type,
    Teuchos::Array<Entity>& adjacent_entities ) const
{
    moab::EntityHandle moab_entity = MoabHelpers::extractEntity(entity);

    std::vector<moab::EntityHandle> adjacencies;
    DTK_CHECK_ERROR_CODE(
	d_moab_mesh->get_moab()->get_adjacencies( 
	    &moab_entity,
	    1,
	    Teuchos::as<int>( entity_type ),
	    true,
	    adjacencies )
	);

    if ( entity.entityType() == entity_type )
    {
	auto remove_it = std::remove( adjacencies.begin(),
				      adjacencies.end(),
				      moab_entity );
	adjacencies.resize( std::distance(adjacencies.begin(),remove_it) );
    }

    int num_adjacencies = adjacencies.size();
    adjacent_entities.resize( num_adjacencies );
    for ( int i = 0; i < num_adjacencies; ++i )
    {
	adjacent_entities[i] = 
	    MoabEntity( adjacencies[i], d_moab_mesh.ptr(), d_set_indexer.ptr() );
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MoabEntitySet.cpp
//---------------------------------------------------------------------------//
