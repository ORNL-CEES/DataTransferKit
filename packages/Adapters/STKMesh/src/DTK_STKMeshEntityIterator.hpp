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
 * \brief DTK_STKMeshEntityIterator.hpp
 * \author Stuart R. Slattery
 * \brief Entity iterator interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STKMESHENTITYITERATOR_HPP
#define DTK_STKMESHENTITYITERATOR_HPP

#include <vector>
#include <functional>

#include "DTK_EntityIterator.hpp"
#include "DTK_Entity.hpp"
#include "DTK_STKMeshEntityIteratorRange.hpp"

#include <Teuchos_RCP.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class STKMeshEntityIterator
  \brief STK mesh entity iterator implementation.
*/
//---------------------------------------------------------------------------//
class STKMeshEntityIterator : public EntityIterator
{
  public:

    /*!
     * \brief Default constructor.
     */
    STKMeshEntityIterator();

    /*! 
     * \brief Constructor.
     */
    STKMeshEntityIterator( 
	const Teuchos::RCP<STKMeshEntityIteratorRange>& entity_range,
	const Teuchos::RCP<stk::mesh::BulkData>& bulk_data,
	const std::function<bool(Entity)>& predicate );

    /*!
     * \brief Copy constructor.
     */
    STKMeshEntityIterator( const STKMeshEntityIterator& rhs );

    /*!
     * \brief Assignment operator.
     */
    STKMeshEntityIterator& operator=( const STKMeshEntityIterator& rhs );

    /*!
     * \brief Destructor.
     */
    ~STKMeshEntityIterator();

    // Pre-increment operator.
    EntityIterator& operator++();

    // Dereference operator.
    Entity& operator*(void);

    // Dereference operator.
    Entity* operator->(void);

    // Equal comparison operator.
    bool operator==( const EntityIterator& rhs ) const;

    // Not equal comparison operator.
    bool operator!=( const EntityIterator& rhs ) const;

    // An iterator assigned to the first valid element in the iterator.
    EntityIterator begin() const;

    // An iterator assigned to the end of all elements under the iterator.
    EntityIterator end() const;

    // Create a clone of the iterator. We need this for the copy constructor
    // and assignment operator to pass along the underlying implementation.
    EntityIterator* clone() const;

  private:

    // Set the current entity.
    void setCurrentEntity();

  private:

    // Range of entities over which the iterator is defined.
    Teuchos::RCP<STKMeshEntityIteratorRange> d_entity_range;

    // Iterator over the entities.
    std::vector<stk::mesh::Entity>::const_iterator d_stk_entity_it;

    // The bulk data owning the entities.
    Teuchos::RCP<stk::mesh::BulkData> d_bulk_data;

    // Current entity.
    Entity d_current_entity;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_STKMESHENTITYITERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntityIterator.hpp
//---------------------------------------------------------------------------//
