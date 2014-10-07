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
 * \brief DTK_AbstractIterator_impl.hpp
 * \author Stuart R. Slattery
 * \brief Abstract iterator interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ABSTRACTITERATOR_IMPL_HPP
#define DTK_ABSTRACTITERATOR_IMPL_HPP

#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class ValueType>
AbstractIterator<ValueType>::AbstractIterator()
{
    // Default predicate always returns true.
    b_predicate = [](ValueType v){return true;};
    b_iterator_impl = NULL;
}

//---------------------------------------------------------------------------//
// Copy constructor.
template<class ValueType>
AbstractIterator<ValueType>::AbstractIterator( 
    const AbstractIterator<ValueType>& rhs )
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
    b_predicate = rhs.b_predicate;
    advanceToFirstValidElement();
}

//---------------------------------------------------------------------------//
// Assignment operator.
template<class ValueType>
AbstractIterator<ValueType>& AbstractIterator<ValueType>::operator=( 
    const AbstractIterator<ValueType>& rhs )
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
    b_predicate = rhs.b_predicate;
    advanceToFirstValidElement();
    return *this;
}

//---------------------------------------------------------------------------//
template<class ValueType>
AbstractIterator<ValueType>::~AbstractIterator()
{
    if ( NULL != b_iterator_impl )
    {
	delete b_iterator_impl;
	b_iterator_impl = NULL;
    }
}

//---------------------------------------------------------------------------//
// Pre-increment operator.
template<class ValueType>
AbstractIterator<ValueType>& AbstractIterator<ValueType>::operator++()
{
    DTK_REQUIRE( NULL != b_iterator_impl );
    DTK_REQUIRE( *b_iterator_impl != b_iterator_impl->end() );
    AbstractIterator<ValueType>& it = b_iterator_impl->operator++();
    AbstractIterator<ValueType> end = b_iterator_impl->end();
    while ( it != end && !b_predicate(*it) )
    {
	it = b_iterator_impl->operator++();
    }
    return it; 
}

//---------------------------------------------------------------------------//
// Post-increment operator.
template<class ValueType>
AbstractIterator<ValueType> AbstractIterator<ValueType>::operator++(ValueType n)
{
    DTK_REQUIRE( NULL != b_iterator_impl );
    DTK_REQUIRE( *b_iterator_impl != b_iterator_impl->end() );
    const AbstractIterator<ValueType> tmp(*this);
    AbstractIterator<ValueType> it = b_iterator_impl->operator++();
    while ( it != b_iterator_impl->end() && !b_predicate(*it) )
    {
	it = b_iterator_impl->operator++();
    }
    return tmp;
}

//---------------------------------------------------------------------------//
// Dereference operator.
template<class ValueType>
ValueType& AbstractIterator<ValueType>::operator*(void)
{
    DTK_REQUIRE( NULL != b_iterator_impl );
    return b_iterator_impl->operator*();
}

//---------------------------------------------------------------------------//
// Dereference operator.
template<class ValueType>
ValueType* AbstractIterator<ValueType>::operator->(void)
{
    DTK_REQUIRE( NULL != b_iterator_impl );
    return b_iterator_impl->operator->();
}

//---------------------------------------------------------------------------//
// Equal comparison operator.
template<class ValueType>
bool AbstractIterator<ValueType>::operator==( 
    const AbstractIterator<ValueType>& rhs ) const
{
    DTK_REQUIRE( NULL != b_iterator_impl );
    return b_iterator_impl->operator==( rhs );
}

//---------------------------------------------------------------------------//
// Not equal comparison operator.
template<class ValueType>
bool AbstractIterator<ValueType>::operator!=( 
    const AbstractIterator<ValueType>& rhs ) const
{
    DTK_REQUIRE( NULL != b_iterator_impl );
    return b_iterator_impl->operator!=( rhs );
}

//---------------------------------------------------------------------------//
// Size of the iterator.
template<class ValueType>
std::size_t AbstractIterator<ValueType>::size() const
{
    if ( NULL != b_iterator_impl )
    {
	return b_iterator_impl->size();
    }
    return 0;
}

//---------------------------------------------------------------------------//
// An iterator assigned to the beginning.
template<class ValueType>
AbstractIterator<ValueType> AbstractIterator<ValueType>::begin() const
{
    if ( NULL != b_iterator_impl )
    {
	return b_iterator_impl->begin();
    }
    return AbstractIterator<ValueType>();
}

//---------------------------------------------------------------------------//
// Create a clone of the iterator.
template<class ValueType>
AbstractIterator<ValueType>* AbstractIterator<ValueType>::clone() const
{
    return new AbstractIterator<ValueType>();
}

//---------------------------------------------------------------------------//
// An iterator assigned to the end.
template<class ValueType>
AbstractIterator<ValueType> AbstractIterator<ValueType>::end() const
{
    if ( NULL != b_iterator_impl )
    {
	return b_iterator_impl->end();
    }
    return AbstractIterator<ValueType>();
}

//---------------------------------------------------------------------------//
// Advance the iterator to the first valid element that satisfies the
// predicate or the end of the iterator. 
template<class ValueType>
void AbstractIterator<ValueType>::advanceToFirstValidElement()
{
    DTK_REQUIRE( NULL != b_iterator_impl );
    if ( !b_predicate(**this) && (*this != end()) )
    {
	operator++();
    }
}

//---------------------------------------------------------------------------//
// Static Members.
//---------------------------------------------------------------------------//
template<class ValueType>
AbstractIterator<ValueType> 
AbstractIterator<ValueType>::setOperation(
    const IteratorSetOperation operation,
    const AbstractIterator<ValueType>& it_1,
    const AbstractIterator<ValueType>& it_2 )
{
    AbstractIterator<ValueType> set_it( it_1 );
    if ( INTERSECTION == operation )
    {
	set_it.b_predicate = std::bind( 
	    std::logical_and<bool>(),
	    std::bind(it_1.b_predicate,std::placeholders::_1),
	    std::bind(it_2.b_predicate,std::placeholders::_1) 
	    );
    }
    else if ( UNION == operation )
    {
    	set_it.b_predicate = std::bind( 
    	    std::logical_or<bool>(),
	    std::bind(it_1.b_predicate,std::placeholders::_1),
	    std::bind(it_2.b_predicate,std::placeholders::_1) 
	    );
    }
    else if ( SUBTRACTION == operation )
    {
    	set_it.b_predicate = std::bind( 
    	    std::logical_and<bool>(),
	    std::bind(it_1.b_predicate,std::placeholders::_1),
	    std::bind(std::logical_not<bool>(),
		      std::bind(it_2.b_predicate,std::placeholders::_1))
	    );
    }
    else
    {
	bool bad_operation_type = true;
	DTK_INSIST( !bad_operation_type );
    }
    set_it.advanceToFirstValidElement();
    return set_it;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ABSTRACTITERATOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_AbstractIterator_impl.hpp
//---------------------------------------------------------------------------//
