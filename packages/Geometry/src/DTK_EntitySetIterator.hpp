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
 * \brief DTK_EntityIterator.hpp
 * \author Stuart R. Slattery
 * \brief Geometric entity set interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYITERATOR_HPP
#define DTK_ENTITYITERATOR_HPP

#include <iterator>
#include <functional>

#include "DTK_GeometricEntity.hpp"

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class EntityIterator
  \brief Geometric entity set iterator interface definition.

  This class provides a mechanism to iterate over a group of entities with a
  specified predicate operation for selection.
*/
//---------------------------------------------------------------------------//
class EntityIterator : 
	public std::iterator<std::forward_iterator_tag,GeometricEntity>
{
  public:

    /*!
     * \brief Constructor.
     */
    EntityIterator();

    /*!
     * \brief Copy constructor.
     */
    EntityIterator( const EntityIterator& rhs );

    /*!
     * \brief Assignment operator.
     */
    EntityIterator& operator=( const EntityIterator& rhs );

    /*!
     * \brief Destructor.
     */
    virtual ~EntityIterator();

    // Pre-increment operator.
    virtual EntityIterator& operator++();

    // Post-increment operator.
    virtual EntityIterator operator++(int);

    // Dereference operator.
    virtual GeometricEntity& operator*(void);

    // Dereference operator.
    virtual GeometricEntity* operator->(void);

    // Equal comparison operator.
    virtual bool operator==( const EntityIterator& rhs ) const;

    // Not equal comparison operator.
    virtual bool operator!=( const EntityIterator& rhs ) const;

    // Create a clone of the iterator. We need this for the copy constructor
    // and assignment operator to pass along the underlying implementation.
    virtual Teuchos::RCP<EntityIterator> clone() const;

    // Size of the iterator.
    virtual std::size_t size() const;

    // An iterator assigned to the beginning.
    virtual EntityIterator begin();

    // An iterator assigned to the end.
    virtual EntityIterator end();

  protected:

    // Implementation.
    Teuchos::RCP<EntityIterator> b_iterator_impl;

    // Predicate.
    std::function<bool(GeometricEntity)> b_predicate;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ENTITYITERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_EntityIterator.hpp
//---------------------------------------------------------------------------//
