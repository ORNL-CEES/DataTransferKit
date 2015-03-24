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
 * \brief DTK_STKMeshEntityIterator.cpp
 * \author Stuart R. Slattery
 * \brief Entity iterator interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_DBC.hpp"
#include "DTK_STKMeshEntityIterator.hpp"
#include <DTK_STKMeshEntity.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Default constructor.
STKMeshEntityIterator::STKMeshEntityIterator()
{
    this->b_iterator_impl = NULL;
}

//---------------------------------------------------------------------------//
// Constructor.
STKMeshEntityIterator::STKMeshEntityIterator(
    const Teuchos::RCP<STKMeshEntityIteratorRange>& entity_range,
    const Teuchos::Ptr<stk::mesh::BulkData>& bulk_data,
    const PredicateFunction& predicate )
    : d_entity_range( entity_range )
    , d_stk_entity_it( d_entity_range->d_stk_entities.begin() )
    , d_bulk_data( bulk_data )
{
    this->b_iterator_impl = NULL;
    this->b_predicate = predicate;
}

//---------------------------------------------------------------------------//
// Copy constructor.
STKMeshEntityIterator::STKMeshEntityIterator( 
    const STKMeshEntityIterator& rhs )
    : d_entity_range( rhs.d_entity_range )
    , d_stk_entity_it( rhs.d_stk_entity_it )
    , d_bulk_data( rhs.d_bulk_data )
{
    this->b_iterator_impl = NULL;
    this->b_predicate = rhs.b_predicate;
}

//---------------------------------------------------------------------------//
// Assignment operator.
STKMeshEntityIterator& STKMeshEntityIterator::operator=( 
    const STKMeshEntityIterator& rhs )
{
    this->b_iterator_impl = NULL;
    this->b_predicate = rhs.b_predicate;
    if ( &rhs == this )
    {
	return *this;
    }
    d_entity_range = rhs.d_entity_range;
    d_stk_entity_it = rhs.d_stk_entity_it;
    d_bulk_data = rhs.d_bulk_data;
    return *this;
}

//---------------------------------------------------------------------------//
// Destructor.
STKMeshEntityIterator::~STKMeshEntityIterator()
{
    this->b_iterator_impl = NULL;
}

//---------------------------------------------------------------------------//
// Pre-increment operator.
EntityIterator& STKMeshEntityIterator::operator++()
{
    ++d_stk_entity_it;
    return *this;
}

//---------------------------------------------------------------------------//
// Dereference operator.
Entity& STKMeshEntityIterator::operator*(void)
{
    this->operator->();
    return d_current_entity;
}

//---------------------------------------------------------------------------//
// Dereference operator.
Entity* STKMeshEntityIterator::operator->(void)
{
    d_current_entity = 
	STKMeshEntity( *d_stk_entity_it, d_bulk_data );
    return &d_current_entity;
}

//---------------------------------------------------------------------------//
// Equal comparison operator.
bool STKMeshEntityIterator::operator==( 
    const EntityIterator& rhs ) const
{ 
    const STKMeshEntityIterator* rhs_it = 
	static_cast<const STKMeshEntityIterator*>(&rhs);
    const STKMeshEntityIterator* rhs_it_impl = 
	static_cast<const STKMeshEntityIterator*>(rhs_it->b_iterator_impl);
    return ( rhs_it_impl->d_stk_entity_it == d_stk_entity_it );
}

//---------------------------------------------------------------------------//
// Not equal comparison operator.
bool STKMeshEntityIterator::operator!=( 
    const EntityIterator& rhs ) const
{
    const STKMeshEntityIterator* rhs_it = 
	static_cast<const STKMeshEntityIterator*>(&rhs);
    const STKMeshEntityIterator* rhs_it_impl = 
	static_cast<const STKMeshEntityIterator*>(rhs_it->b_iterator_impl);
    return ( rhs_it_impl->d_stk_entity_it != d_stk_entity_it );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the beginning.
EntityIterator STKMeshEntityIterator::begin() const
{ 
    return STKMeshEntityIterator( 
	d_entity_range, d_bulk_data, this->b_predicate );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the end.
EntityIterator STKMeshEntityIterator::end() const
{
    STKMeshEntityIterator end_it( 
	d_entity_range, d_bulk_data, this->b_predicate );
    end_it.d_stk_entity_it = d_entity_range->d_stk_entities.end();
    return end_it;
}

//---------------------------------------------------------------------------//
// Create a clone of the iterator. We need this for the copy constructor
// and assignment operator to pass along the underlying implementation.
EntityIterator* STKMeshEntityIterator::clone() const
{
    return new STKMeshEntityIterator(*this);
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntityIterator.cpp
//---------------------------------------------------------------------------//
