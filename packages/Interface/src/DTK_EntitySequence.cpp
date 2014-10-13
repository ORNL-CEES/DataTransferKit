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
 * \brief DTK_EntitySequence.cpp
 * \author Stuart R. Slattery
 * \brief Entity sequence interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_EntitySequence.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
EntitySequence::EntitySequence()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
EntitySequence::~EntitySequence()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Number of elements in the sequence that meet the predicate criteria.
std::size_t size() const
{
    return std::distance( this->begin(), this->end() );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the first valid element in the sequence.
ForwardIterator EntitySequence::begin() const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return ForwardIterator();
}

//---------------------------------------------------------------------------//
// An iterator assigned to the end of all elements under the iterator.
ForwardIterator EntitySequence::end() const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return ForwardIterator();
}

//---------------------------------------------------------------------------//
// A filtered iterator assigned to the first valid element in the sequence.
FilteredForwardIterator 
EntitySequence::filteredBegin( const Predicate& predicate ) const
{
    return FilteredForwardIterator( predicate, this->begin(), this->end() );
}

//---------------------------------------------------------------------------//
// An filtered iterator  assigned to the last valid element in the sequence.
FilteredForwardIterator 
EntitySequence::filteredEnd( const Predicate& predicate ) const
{
    return FilteredForwardIterator( predicate, this->end(), this->end() );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_EntitySequence.cpp
//---------------------------------------------------------------------------//
