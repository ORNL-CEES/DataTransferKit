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
 * \brief DTK_EntityIterator.cpp
 * \author Stuart R. Slattery
 * \brief Entity iterator interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_DBC.hpp"
#include "DTK_EntityIterator.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
EntityIterator::EntityIterator()
{
    b_iterator_impl = NULL;
}

//---------------------------------------------------------------------------//
// Copy constructor.
EntityIterator::EntityIterator( 
    const EntityIterator& rhs )
{
    b_iterator_impl = NULL;
    if ( NULL == rhs.b_iterator_impl )
    {
	b_iterator_impl = rhs.clone();
    }
    else
    {
	b_iterator_impl = rhs.b_iterator_impl->clone();
    }
}

//---------------------------------------------------------------------------//
// Assignment operator.
EntityIterator& EntityIterator::operator=( 
    const EntityIterator& rhs )
{
    if ( this == &rhs )
    {
	return *this;
    }
    if ( NULL != b_iterator_impl )
    {
	delete b_iterator_impl;
	b_iterator_impl = NULL;
    }
    if ( NULL == rhs.b_iterator_impl )
    {
	b_iterator_impl = rhs.clone();
    }
    else
    {
	b_iterator_impl = rhs.b_iterator_impl->clone();
    }
    return *this;
}

//---------------------------------------------------------------------------//
EntityIterator::~EntityIterator()
{
    if ( NULL != b_iterator_impl )
    {
	delete b_iterator_impl;
	b_iterator_impl = NULL;
    }
}

//---------------------------------------------------------------------------//
// Pre-increment operator.
EntityIterator& EntityIterator::operator++()
{
    DTK_REQUIRE( NULL != b_iterator_impl );
    DTK_REQUIRE( *b_iterator_impl != b_iterator_impl->end() );
    return b_iterator_impl->operator++();
}

//---------------------------------------------------------------------------//
// Post-increment operator.
EntityIterator EntityIterator::operator++(int n)
{
    DTK_REQUIRE( NULL != b_iterator_impl );
    DTK_REQUIRE( *b_iterator_impl != b_iterator_impl->end() );
    return b_iterator_impl->operator++(n);
}

//---------------------------------------------------------------------------//
// Dereference operator.
Entity& EntityIterator::operator*(void)
{
    DTK_REQUIRE( NULL != b_iterator_impl );
    return b_iterator_impl->operator*();
}

//---------------------------------------------------------------------------//
// Dereference operator.
Entity* EntityIterator::operator->(void)
{
    DTK_REQUIRE( NULL != b_iterator_impl );
    return b_iterator_impl->operator->();
}

//---------------------------------------------------------------------------//
// Equal comparison operator.
bool EntityIterator::operator==( 
    const EntityIterator& rhs ) const
{
    DTK_REQUIRE( NULL != b_iterator_impl );
    return b_iterator_impl->operator==( rhs );
}

//---------------------------------------------------------------------------//
// Not equal comparison operator.
bool EntityIterator::operator!=( 
    const EntityIterator& rhs ) const
{
    DTK_REQUIRE( NULL != b_iterator_impl );
    return b_iterator_impl->operator!=( rhs );
}

//---------------------------------------------------------------------------//
// Size of the iterator. This is the number of objects in the iterator that
// meet the predicate criteria.
std::size_t EntityIterator::size() const
{
    if ( NULL != b_iterator_impl )
    {
	return b_iterator_impl->size();
    }
    return 0;
}

//---------------------------------------------------------------------------//
// An iterator assigned to the beginning.
EntityIterator EntityIterator::begin() const
{
    return b_iterator_impl->begin();
}

//---------------------------------------------------------------------------//
// An iterator assigned to the end.
EntityIterator EntityIterator::end() const
{
   return b_iterator_impl->end();
}

//---------------------------------------------------------------------------//
// Create a clone of the iterator.
EntityIterator* EntityIterator::clone() const
{
    return new EntityIterator();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_EntityIterator.cpp
//---------------------------------------------------------------------------//
