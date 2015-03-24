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
 * \brief DTK_BasicEntitySet.hpp
 * \author Stuart R. Slattery
 * \brief Basic entity set implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BASICENTITYSET_HPP
#define DTK_BASICENTITYSET_HPP

#include <unordered_map>

#include "DTK_Types.hpp"
#include "DTK_EntitySet.hpp"
#include "DTK_Entity.hpp"
#include "DTK_EntityIterator.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class BasicEntitySetIterator
  \brief ementation of iterator over entities in a basic set.
*/
class BasicEntitySetIterator : public EntityIterator
{
  public:

    // Default constructor.
    BasicEntitySetIterator();

    // Constructor.
    BasicEntitySetIterator( 
	Teuchos::RCP<std::unordered_map<EntityId,Entity> > map,
	const PredicateFunction& predicate );

    // Copy constructor.
    BasicEntitySetIterator( const BasicEntitySetIterator& rhs );

    // Destructor.
    ~BasicEntitySetIterator();
    
    /*!
     * \brief Assignment operator.
     */
    BasicEntitySetIterator& operator=( const BasicEntitySetIterator& rhs );

    // Pre-increment operator.
    EntityIterator& operator++() override;

    // Dereference operator.
    Entity& operator*(void) override;

    // Dereference operator.
    Entity* operator->(void) override;

    // Equal comparison operator.
    bool operator==( const EntityIterator& rhs ) const override;

    // Not equal comparison operator.
    bool operator!=( const EntityIterator& rhs ) const override;

    // An iterator assigned to the beginning.
    EntityIterator begin() const override;

    // An iterator assigned to the end.
    EntityIterator end() const override;

    // Create a clone of the iterator. We need this for the copy constructor
    // and assignment operator to pass along the underlying implementation.
    EntityIterator* clone() const override;

  private:

    // Map to iterate over.
    Teuchos::RCP<std::unordered_map<EntityId,Entity> > d_map;

    // Iterator over the entity map.
    std::unordered_map<EntityId,Entity>::iterator d_map_it;

    // Pointer to the current entity.
    Entity* d_entity;
};

//---------------------------------------------------------------------------//
/*!
  \class BasicEntitySet
  \brief Basic implementation of the entity set interface.
*/
//---------------------------------------------------------------------------//
class BasicEntitySet : public EntitySet
{
  public:

    //@{
    //! Entity map typedefs.
    typedef std::pair<EntityId,Entity> EntityIdPair;
    //@}

    /*!
     * \brief Constructor.
     */
    BasicEntitySet( const Teuchos::RCP<const Teuchos::Comm<int> > comm,
		    const int physical_dimension );

    /*!
     * \brief Add an entity to the set.
     */
    void addEntity( const Entity& entity );

    //@{
    //! Parallel functions.
    /*!
     * \brief Get the parallel communicator for the entity set.
     * \return A reference-counted pointer to the parallel communicator.
     */
    Teuchos::RCP<const Teuchos::Comm<int> > communicator() const override;
    //@}

    //@{
    //! Geometric data functions.
    /*!
     * \brief Return the physical dimension of the entities in the set.
     * \return The physical dimension of the set.
     */
    int physicalDimension() const override;
    //@}

    //@{
    //! Entity access functions.
    /*!
     * \brief Given an EntityId, get the entity.
     * \param entity_id Get the entity with this id.
     * \param entity The entity with the given id.
     */
    void getEntity( const EntityType entity_type,
		    const EntityId entity_id, 
		    Entity& entity ) const override;

    /*!
     * \brief Get an iterator over a subset of the entity set that satisfies
     * the given predicate.
     * \param entity_type The type of entity to get an iterator over.
     * \param predicate A predicate to select the entity set to iterate over.
     * \return An iterator to the entities that satisfy the predicate.
     */
    EntityIterator entityIterator(
	const EntityType entity_type,
	const PredicateFunction& predicate ) const override;

    /*!
     * \brief Given an entity, get the entities of the given type that are
     * adjacent to it.
     */
    virtual void getAdjacentEntities(
	const Entity& entity,
	const EntityType entity_type,
	Teuchos::Array<Entity>& adjacent_entities ) const override;
    //@}

  private:

    // Given a parametric dimension, get an id-to-entity map.
    std::unordered_map<EntityId,Entity>& getEntityMap(
	const int parametric_dimension );
    
  private:

    // Parallel communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > d_comm;

    // Physical dimension.
    int d_physical_dim;

    // Id-to-entity maps.
    mutable Teuchos::Array<std::unordered_map<EntityId,Entity> > d_entities;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_BASICENTITYSET_HPP

//---------------------------------------------------------------------------//
// end DTK_BasicEntitySet.hpp
//---------------------------------------------------------------------------//
