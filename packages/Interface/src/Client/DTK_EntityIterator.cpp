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
 * \brief DTK_EntityIterator.cpp
 * \author Stuart R. Slattery
 * \brief Entity iterator interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_EntityIterator.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
EntityIterator::EntityIterator()
    : b_iterator_impl( nullptr )
{
    // Default predicate always returns true.
    b_predicate = []( Entity v ) { return true; };
}

//---------------------------------------------------------------------------//
// Copy constructor.
EntityIterator::EntityIterator( const EntityIterator &rhs )
{
    b_iterator_impl.reset();
    if ( rhs.b_iterator_impl )
    {
        b_iterator_impl = std::move( rhs.b_iterator_impl->clone() );
        b_predicate = b_iterator_impl->b_predicate;
    }
    else
    {
        b_iterator_impl = std::move( rhs.clone() );
        b_predicate = rhs.b_predicate;
    }
}

//---------------------------------------------------------------------------//
// Assignment operator.
EntityIterator &EntityIterator::operator=( const EntityIterator &rhs )
{
    b_iterator_impl.reset();
    if ( this == &rhs )
    {
        return *this;
    }
    if ( rhs.b_iterator_impl )
    {
        b_iterator_impl = std::move( rhs.b_iterator_impl->clone() );
        b_predicate = b_iterator_impl->b_predicate;
    }
    else
    {
        b_iterator_impl = std::move( rhs.clone() );
        b_predicate = rhs.b_predicate;
    }
    return *this;
}

//---------------------------------------------------------------------------//
EntityIterator::~EntityIterator() { /* ... */}

//---------------------------------------------------------------------------//
// Pre-increment operator.
EntityIterator &EntityIterator::operator++()
{
    DTK_REQUIRE( b_iterator_impl );
    DTK_REQUIRE( *b_iterator_impl != b_iterator_impl->end() );

    increment();
    return *b_iterator_impl;
}

//---------------------------------------------------------------------------//
// Post-increment operator.
EntityIterator EntityIterator::operator++( int n )
{
    DTK_REQUIRE( b_iterator_impl );
    DTK_REQUIRE( *b_iterator_impl != b_iterator_impl->end() );

    const EntityIterator tmp( *this );
    increment();
    return tmp;
}

//---------------------------------------------------------------------------//
// Dereference operator.
Entity &EntityIterator::operator*( void )
{
    DTK_REQUIRE( b_iterator_impl );
    return b_iterator_impl->operator*();
}

//---------------------------------------------------------------------------//
// Dereference operator.
Entity *EntityIterator::operator->( void )
{
    DTK_REQUIRE( b_iterator_impl );
    return b_iterator_impl->operator->();
}

//---------------------------------------------------------------------------//
// Equal comparison operator.
bool EntityIterator::operator==( const EntityIterator &rhs ) const
{
    if ( nullptr == b_iterator_impl )
        return ( nullptr == b_iterator_impl );
    return b_iterator_impl->operator==( rhs );
}

//---------------------------------------------------------------------------//
// Not equal comparison operator.
bool EntityIterator::operator!=( const EntityIterator &rhs ) const
{
    return !( b_iterator_impl->operator==( rhs ) );
}

//---------------------------------------------------------------------------//
// Size of the iterator. This is the number of objects in the iterator that
// meet the predicate criteria.
std::size_t EntityIterator::size() const
{
    std::size_t size = 0;
    if ( b_iterator_impl )
    {
        size = std::distance( this->begin(), this->end() );
    }
    return size;
}

//---------------------------------------------------------------------------//
// An iterator assigned to the beginning.
EntityIterator EntityIterator::begin() const
{
    EntityIterator begin_it;
    if ( b_iterator_impl )
    {
        begin_it = b_iterator_impl->begin();
        begin_it.b_predicate = b_predicate;
        begin_it.advanceToFirstValidElement();
    }
    return begin_it;
}

//---------------------------------------------------------------------------//
// An iterator assigned to the end.
EntityIterator EntityIterator::end() const
{
    if ( b_iterator_impl )
    {
        return b_iterator_impl->end();
    }
    return EntityIterator();
}

//---------------------------------------------------------------------------//
// Create a clone of the iterator.
std::unique_ptr<EntityIterator> EntityIterator::clone() const
{
    return std::unique_ptr<EntityIterator>( new EntityIterator() );
}

//---------------------------------------------------------------------------//
// Advance the iterator to the first valid element that satisfies the
// predicate or the end of the iterator.
void EntityIterator::advanceToFirstValidElement()
{
    DTK_REQUIRE( b_iterator_impl );
    if ( ( *this != end() ) && !b_predicate( **this ) )
    {
        increment();
    }
}

//---------------------------------------------------------------------------//
// Increment the iterator implementation forward until either a valid
// increment is found or we have reached the end.
void EntityIterator::increment()
{
    DTK_REQUIRE( b_iterator_impl );
    DTK_REQUIRE( *b_iterator_impl != b_iterator_impl->end() );

    // Apply the increment operator.
    EntityIterator &it = b_iterator_impl->operator++();

    // Get the end of the range.
    EntityIterator end = b_iterator_impl->end();

    // If the we are not at the end or the predicate is not satisfied by the
    // current element, increment until either of these conditions is
    // satisfied.
    while ( it != end && !b_predicate( *it ) )
    {
        it = b_iterator_impl->operator++();
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_EntityIterator.cpp
//---------------------------------------------------------------------------//
