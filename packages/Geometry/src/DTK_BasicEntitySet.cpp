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
 * \brief DTK_BasicEntitySet.cpp
 * \author Stuart R. Slattery
 * \brief Basic entity set implementation.
 */
//---------------------------------------------------------------------------//

#include "DTK_BasicEntitySet.hpp"
#include "DTK_DBC.hpp"
#include "DTK_Box.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_OrdinalTraits.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Default constructor.
BasicEntitySet::BasicEntitySet()
    : d_entities( 4 )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Constructor.
BasicEntitySet::BasicEntitySet(
    const Teuchos::RCP<const Teuchos::Comm<int> > comm,
    const int physical_dimension )
    : d_comm( comm )
    , d_dimension( physical_dimension )
    , d_entities( 4 )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
BasicEntitySet::~BasicEntitySet()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Return a string indicating the derived entity set type.
std::string BasicEntitySet::entitySetType() const
{
    return std::string("DTK Basic Entity Set");
}

//---------------------------------------------------------------------------//
// Assign a parallel communicator to the entity set. This will only be done
// immediately after construct through the AbstractBuilder interface.
void BasicEntitySet::assignCommunicator( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
{
    d_comm = comm;
}

//---------------------------------------------------------------------------//
// Get the parallel communicator for the entity set.
Teuchos::RCP<const Teuchos::Comm<int> >
BasicEntitySet::communicator() const
{
    return d_comm;
}

//---------------------------------------------------------------------------//
// Get the local number of entities in the set of the given
std::size_t BasicEntitySet::localNumberOfEntities( 
    const int parametric_dimension ) const
{
    DTK_REQUIRE( parametric_dimension <= d_dimension );
    return d_entities[parametric_dimension].size();
}

//---------------------------------------------------------------------------//
// Get the global number of entities in the set of the given
std::size_t BasicEntitySet::globalNumberOfEntities(
    const int parametric_dimension ) const
{
    DTK_REQUIRE( parametric_dimension <= d_dimension );
    std::size_t num_local = localNumberOfEntities( parametric_dimension );
    std::size_t num_global = 0;
    Teuchos::reduceAll( *d_comm, Teuchos::REDUCE_SUM, 
			num_local, Teuchos::Ptr<std::size_t>(&num_global) );
    return num_global;
}
    
//---------------------------------------------------------------------------//
// Get a forward iterator assigned to the beginning of the entities in
// the set of the given parametric dimension. 
std::iterator<std::forward_iterator_tag,GeometricEntity>
BasicEntitySet::entityIteratorBegin( const int parametric_dimension ) const
{
    DTK_REQUIRE( parametric_dimension <= d_dimension );
    std::unordered_map<EntityId,GeometricEntity>::const_iterator begin =
	d_entities[parametric_dimension].begin();
    return BasicEntitySetIterator( begin );
}

//---------------------------------------------------------------------------//
// Get a forward iterator assigned to the end of the entities in the set
// of the given parametric dimension.
std::iterator<std::forward_iterator_tag,GeometricEntity>
BasicEntitySet::entityIteratorEnd( const int parametric_dimension ) const
{
    DTK_REQUIRE( parametric_dimension <= d_dimension );
    std::unordered_map<EntityId,GeometricEntity>::const_iterator end =
	d_entities[parametric_dimension].end();
    return BasicEntitySetIterator( end );
}

//---------------------------------------------------------------------------//
// Get the identifiers for all local entities in the set of a given
void BasicEntitySet::localEntityIds( 
    const int parametric_dimension,
    const Teuchos::ArrayView<EntityId>& entity_ids ) const
{
    DTK_REQUIRE( parametric_dimension <= d_dimension );
    DTK_REQUIRE( Teuchos::as<std::size_t>(entity_ids.size()) == 
		 localNumberOfEntities(parametric_dimension) );
    Teuchos::ArrayView<EntityId>::iterator id_it = entity_ids.begin();

    std::iterator<std::forward_iterator_tag,GeometricEntity> entity_begin = 
	entityIteratorBegin( parametric_dimension );
    std::iterator<std::forward_iterator_tag,GeometricEntity> entity_end =
	entityIteratorEnd( parametric_dimension );

    std::iterator<std::forward_iterator_tag,GeometricEntity> entity_it;
    for ( entity_it = entity_begin;
	  entity_it != entity_end;
	  ++entity_it )
    {
	*id_it = entity_it->id();
	++id_it;
    }
}

//---------------------------------------------------------------------------//
// Given an EntityId, get the entity.
void BasicEntitySet::getEntity( 
    const EntityId entity_id, GeometricEntity& entity ) const
{
    DTK_REQUIRE( d_entity_dims.count(entity_id) );
    int entity_dim = d_entity_dims.find(entity_id)->second;
    DTK_CHECK( d_entities[entity_dim].count(entity_id) );
    entity = d_entities[entity_dim].find(entity_id)->second;
}

//---------------------------------------------------------------------------//
// Indicate that the entity set will be modified.
void BasicEntitySet::startModification()
{
    // We can always modify the basic entity set.
}

//---------------------------------------------------------------------------//
// Add an entity to the set.
void BasicEntitySet::addEntity( const GeometricEntity& entity )
{
    int parametric_dimension = entity.parametricDimension();
    DTK_CHECK( parametric_dimension <= d_dimension );
    d_entity_dims.insert( 
	std::pair<EntityId,int>(entity.id(),parametric_dimension) );
    d_entities[parametric_dimension].insert( EntityIdPair(entity.id(),entity) );
}

//---------------------------------------------------------------------------//
// Indicate that modification of the entity set is complete.
void BasicEntitySet::endModification()
{
    /* ... */
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entities in the set.
int BasicEntitySet::physicalDimension() const
{
    return d_dimension;
}

//---------------------------------------------------------------------------//
// Get the local bounding box of entities of the set.
void BasicEntitySet::localBoundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    double max = std::numeric_limits<double>::max();
    Box bounding_box( d_comm->getRank(), d_comm->getRank(),
		      max, max, max, -max, -max, -max );
    Box entity_box;
    std::iterator<std::forward_iterator_tag,GeometricEntity> entity_begin;
    std::iterator<std::forward_iterator_tag,GeometricEntity> entity_end;
    std::iterator<std::forward_iterator_tag,GeometricEntity> entity_it;
    Teuchos::Tuple<double,6> entity_bounds;
    for ( int i = 0; i < 4; ++i )
    {
	entity_begin = entityIteratorBegin( i );
	entity_end = entityIteratorEnd( i );
	for ( entity_it = entity_begin;
	      entity_it != entity_end;
	      ++entity_it )
	{
	    entity_it->boundingBox( entity_bounds );
	    entity_box = Box( 0, 0, entity_bounds );
	    bounding_box += entity_box;
	}
    }
    bounds = bounding_box.getBounds();
}

//---------------------------------------------------------------------------//
// Get the global bounding box of entities of the set.
void BasicEntitySet::globalBoundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    Teuchos::Tuple<double,6> local_bounds;
    localBoundingBox( local_bounds );
    Box local_box = Box( 0, 0, local_bounds );
    Box bounding_box = Box();
    Teuchos::reduceAll( *d_comm, Teuchos::REDUCE_SUM, 
			local_box, Teuchos::Ptr<Box>(&bounding_box) );
    bounds = bounding_box.getBounds();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_BasicEntitySet.cpp
//---------------------------------------------------------------------------//
