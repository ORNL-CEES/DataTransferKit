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
 * \brief DTK_EntitySetIterator.hpp
 * \author Stuart R. Slattery
 * \brief Geometric entity set interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYSETITERATOR_HPP
#define DTK_ENTITYSETITERATOR_HPP

#include <iterator>

#include "DTK_GeometricEntity.hpp"

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class EntitySetIterator
  \brief Geometric entity set iterator interface definition.
*/
//---------------------------------------------------------------------------//
class EntitySetIterator : 
	public std::iterator<std::forward_iterator_tag,GeometricEntity>
{
  public:

    /*!
     * \brief Constructor.
     */
    EntitySetIterator();

    /*!
     * \brief Destructor.
     */
    virtual ~EntitySetIterator();

    // Pre-increment operator.
    virtual EntitySetIterator& operator++();

    // Post-increment operator.
    virtual EntitySetIterator operator++(int);

    // Dereference operator.
    virtual GeometricEntity& operator*(void);

    // Dereference operator.
    virtual GeometricEntity* operator->(void);

    // Equal comparison operator.
    virtual bool operator==( const EntitySetIterator& rhs ) const;

    // Not equal comparison operator.
    virtual bool operator!=( const EntitySetIterator& rhs ) const;

  protected:

    // Implementation.
    Teuchos::RCP<EntitySetIterator> b_iterator_impl;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ENTITYSETITERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_EntitySetIterator.hpp
//---------------------------------------------------------------------------//
