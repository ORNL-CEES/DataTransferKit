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
 * \brief DTK_BasicEntitySetImplementation.cpp
 * \author Stuart R. Slattery
 * \brief Basic entity set implementation.
 */
//---------------------------------------------------------------------------//

#include "DTK_BasicEntitySetImplementation.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Ptr.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
BasicEntitySetImplementation::BasicEntitySetImplementation(
    const Teuchos::RCP<const Teuchos::Comm<int> > comm,
    const int physical_dimension )
    : d_comm( comm )
    , d_dimension( physical_dimension )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
BasicEntitySetImplementation::~BasicEntitySetImplementation()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Return a string indicating the derived entity set type.
std::string BasicEntitySetImplementation::entitySetType() const
{
    return std::string("DTK Basic Entity Set");
}

//---------------------------------------------------------------------------//
// Get the parallel communicator for the entity set.
Teuchos::RCP<const Teuchos::Comm<int> >
BasicEntitySetImplementation::communicator() const
{
    return d_comm;
}

//---------------------------------------------------------------------------//
// Get the local number of entities in the set of the given
std::size_t BasicEntitySetImplementation::localNumberOfEntities( 
    const int parametric_dimension ) const
{
    std::size_t num_local = 0;
    ConstEntityIterator entity_it;
    for ( entity_it = d_entities.begin();
	  entity_it != d_entities.end();
	  ++entity_it )
    {
	if ( entity_it->second->parametricDimension() == parametric_dimension )
	{
	    ++num_local;
	}
    }
    return num_local;
}

//---------------------------------------------------------------------------//
// Get the global number of entities in the set of the given
std::size_t BasicEntitySetImplementation::globalNumberOfEntities(
    const int parametric_dimension ) const
{
    std::size_t num_local = localNumberOfEntities( parameteric_dimension );
    std::size_t num_global = 0;
    Teuchos::reduceAll( *d_comm, Teuchos::REDUCE_SUM, 
			num_local, Teuchos::Ptr<std::size_t>(&num_global) );
    return num_global;
}
    
//---------------------------------------------------------------------------//
// Get the identifiers for all local entities in the set of a given
void BasicEntitySetImplementation::localEntityIds( 
    const int parametric_dimension,
    const Teuchos::ArrayView<EntityId>& entity_ids ) const
{
    DTK_REQUIRE( entity_ids.size() == 
		 localNumberOfEntities(parametric_dimension) );
    Teuchos::ArrayView<EntityId>::const_iterator id_it = entity_ids.begin();
    ConstEntityIterator entity_it;
    for ( entity_it = d_entities.begin();
	  entity_it != d_entities.end();
	  ++entity_it )
    {
	if ( entity_it->second->parametricDimension() == parametric_dimension )
	{
	    *id_it = entity_it->second->id();
	    ++id_it;
	}
    }
}

//---------------------------------------------------------------------------//
// Given an EntityId, get the entity.
void BasicEntitySetImplementation::getEntity( 
    const EntityId entity_id, 
    Teuchos::RCP<GeometricEntity>& entity ) const
{
    DTK_REQUIRE( d_entities.count(entity_id) );
    return d_entities.find(entity_id)->second;
}

//---------------------------------------------------------------------------//
// Indicate that the entity set will be modified.
void BasicEntitySetImplementation::startModification()
{
    // We can always modify the entity set.
}

//---------------------------------------------------------------------------//
// Add an entity to the set.
void BasicEntitySetImplementation::addEntity( 
    const Teuchos::RCP<GeometricEntity>& entity )
{
    DTK_CHECK( Teuchos::nonnull(entity) );
    d_entities.insert( EntityIdPair(entity->id(),entity) );
}

//---------------------------------------------------------------------------//
// Indicate that modification of the entity set is complete.
void BasicEntitySetImplementation::endModification()
{
    /* ... */
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entities in the set.
int BasicEntitySetImplementation::physicalDimension() const
{
    return d_dimension;
}

//---------------------------------------------------------------------------//
// Get the local bounding box of entities of the set.
void BasicEntitySetImplementation::localBoundingBox( Box& bounding_box ) const
{
    bounding_box = Box( comm->getRank(), comm->getRank(),
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
    ConstEntityIterator entity_it;
    for ( entity_it = d_entities.begin();
	  entity_it != d_entities.end();
	  ++entity_it )
    {
	bounding_box += entity_it->second->boundingBox();
    }
}

//---------------------------------------------------------------------------//
// Get the global bounding box of entities of the set.
void BasicEntitySetImplementation::globalBoundingBox( Box& bounding_box ) const
{
    Box local_box;
    localBoundingBox( local_box );
    bounding_box = Box();
    Teuchos::reduceAll( *d_comm, Teuchos::REDUCE_SUM, 
			local_box, Teuchos::Ptr<Box>(&bounding_box) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_BasicEntitySetImplementation.cpp
//---------------------------------------------------------------------------//
