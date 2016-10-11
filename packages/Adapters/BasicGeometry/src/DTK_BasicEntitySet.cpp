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
 * \brief DTK_BasicEntitySet.cpp
 * \author Stuart R. Slattery
 * \brief Basic entity set implementation.
 */
//---------------------------------------------------------------------------//

#include "DTK_BasicEntitySet.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Ptr.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// BasicEntitySetIterator implementation.
//---------------------------------------------------------------------------//
// Default constructor.
BasicEntitySetIterator::BasicEntitySetIterator()
    : d_entity( NULL )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Constructor.
BasicEntitySetIterator::BasicEntitySetIterator(
    Teuchos::RCP<std::unordered_map<EntityId, Entity>> map,
    const PredicateFunction &predicate )
    : d_map( map )
    , d_map_it( d_map->begin() )
{
    if ( ( d_map->size() > 0 ) && ( d_map_it != d_map->end() ) )
    {
        d_entity = &( d_map_it->second );
    }
    this->b_predicate = predicate;
}

//---------------------------------------------------------------------------//
// Copy constructor.
BasicEntitySetIterator::BasicEntitySetIterator(
    const BasicEntitySetIterator &rhs )
    : d_map( rhs.d_map )
    , d_map_it( rhs.d_map_it )
{
    if ( ( d_map->size() > 0 ) && ( d_map_it != d_map->end() ) )
    {
        d_entity = &( d_map_it->second );
    }
    this->b_predicate = rhs.b_predicate;
}

//---------------------------------------------------------------------------//
// Assignment operator.
BasicEntitySetIterator &BasicEntitySetIterator::
operator=( const BasicEntitySetIterator &rhs )
{
    this->b_predicate = rhs.b_predicate;
    if ( &rhs == this )
    {
        return *this;
    }
    d_map = rhs.d_map;
    d_map_it = rhs.d_map_it;
    if ( ( d_map->size() > 0 ) && ( d_map_it != d_map->end() ) )
    {
        d_entity = &( d_map_it->second );
    }
    return *this;
}

//---------------------------------------------------------------------------//
// Pre-increment operator.
EntityIterator &BasicEntitySetIterator::operator++()
{
    ++d_map_it;
    return *this;
}

//---------------------------------------------------------------------------//
// Dereference operator.
Entity &BasicEntitySetIterator::operator*( void )
{
    this->operator->();
    return *d_entity;
}

//---------------------------------------------------------------------------//
// Dereference operator.
Entity *BasicEntitySetIterator::operator->( void )
{
    DTK_REQUIRE( d_map_it != d_map->end() );
    d_entity = &( d_map_it->second );
    return d_entity;
}

//---------------------------------------------------------------------------//
// Equal comparison operator.
bool BasicEntitySetIterator::operator==( const EntityIterator &rhs ) const
{
    const BasicEntitySetIterator *rhs_vec =
        static_cast<const BasicEntitySetIterator *>( &rhs );
    const BasicEntitySetIterator *rhs_vec_impl =
        static_cast<const BasicEntitySetIterator *>(
            rhs_vec->b_iterator_impl.get() );
    return ( rhs_vec_impl->d_map_it == d_map_it );
}

//---------------------------------------------------------------------------//
// Not equal comparison operator.
bool BasicEntitySetIterator::operator!=( const EntityIterator &rhs ) const
{
    const BasicEntitySetIterator *rhs_vec =
        static_cast<const BasicEntitySetIterator *>( &rhs );
    const BasicEntitySetIterator *rhs_vec_impl =
        static_cast<const BasicEntitySetIterator *>(
            rhs_vec->b_iterator_impl.get() );
    return ( rhs_vec_impl->d_map_it != d_map_it );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the beginning.
EntityIterator BasicEntitySetIterator::begin() const
{
    return BasicEntitySetIterator( d_map, this->b_predicate );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the end.
EntityIterator BasicEntitySetIterator::end() const
{
    BasicEntitySetIterator end_it( d_map, this->b_predicate );
    end_it.d_map_it = d_map->end();
    return end_it;
}

//---------------------------------------------------------------------------//
// Create a clone of the iterator. We need this for the copy constructor
// and assignment operator to pass along the underlying implementation.
std::unique_ptr<EntityIterator> BasicEntitySetIterator::clone() const
{
    return std::unique_ptr<EntityIterator>(
        new BasicEntitySetIterator( *this ) );
}

//---------------------------------------------------------------------------//
// BasicEntitySet implementation.
//---------------------------------------------------------------------------//
// Constructor.
BasicEntitySet::BasicEntitySet(
    const Teuchos::RCP<const Teuchos::Comm<int>> comm,
    const int physical_dimension )
    : d_comm( comm )
    , d_physical_dim( physical_dimension )
    , d_entities( 4 )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Add an entity to the set.
void BasicEntitySet::addEntity( const Entity &entity )
{
    d_entities[entity.topologicalDimension()].emplace( entity.id(), entity );
}

//---------------------------------------------------------------------------//
// Get the parallel communicator for the entity set.
Teuchos::RCP<const Teuchos::Comm<int>> BasicEntitySet::communicator() const
{
    return d_comm;
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entities in the set.
int BasicEntitySet::physicalDimension() const { return d_physical_dim; }

//---------------------------------------------------------------------------//
// Given an EntityId, get the entity.
void BasicEntitySet::getEntity( const EntityId entity_id,
                                const int topological_dimension,
                                Entity &entity ) const
{
    DTK_CHECK( d_entities[topological_dimension].count( entity_id ) );
    entity = d_entities[topological_dimension].find( entity_id )->second;
}

//---------------------------------------------------------------------------//
// Get an iterator over a subset of the entity set that satisfies the given
// predicate.
EntityIterator
BasicEntitySet::entityIterator( const int topological_dimension,
                                const PredicateFunction &predicate ) const
{
    Teuchos::RCP<std::unordered_map<EntityId, Entity>> map_ptr =
        Teuchos::rcpFromRef( d_entities[topological_dimension] );
    return BasicEntitySetIterator( map_ptr, predicate );
}

//---------------------------------------------------------------------------//
// Given an entity, get the entities of the given type that are adjacent to
// it.
void BasicEntitySet::getAdjacentEntities(
    const Entity &entity, const int adjacent_dimension,
    Teuchos::Array<Entity> &adjacent_entities ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_BasicEntitySet.cpp
//---------------------------------------------------------------------------//
