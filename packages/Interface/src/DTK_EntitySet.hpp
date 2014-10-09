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
 * \brief DTK_EntitySet.hpp
 * \author Stuart R. Slattery
 * \brief Geometric entity set interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYSET_HPP
#define DTK_ENTITYSET_HPP

#include <functional>

#include "DTK_Types.hpp"
#include "DTK_Entity.hpp"
#include "DTK_AbstractBuilder.hpp"
#include "DTK_AbstractBuildableObject.hpp"
#include "DTK_AbstractIterator.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class EntitySet
  \brief Geometric entity set interface definition.

  An entity set is a coherent collection of geometric entities with a parallel
  distribution.
*/
//---------------------------------------------------------------------------//
class EntitySet : public AbstractBuildableObject<EntitySet>
{
  public:

    /*!
     * \brief Constructor.
     */
    EntitySet();

    /*!
     * \brief Destructor.
     */
    virtual ~EntitySet();

    //@{
    //! Identification functions.
    /*!
     * \brief Return a string indicating the derived entity set type.
     * \return A string key for the derived set type.
     */
    virtual std::string name() const;
    //@}

    //@{
    //! Parallel functions.
    /*!
     * \brief Assign a parallel communicator to the entity set. This will only
     * be done immediately after construct through the AbstractBuilder
     * interface.
     * \param comm The communicator to set with the entity set.
     */
    virtual void assignCommunicator(
	const Teuchos::RCP<const Teuchos::Comm<int> >& comm );

    /*!
     * \brief Get the parallel communicator for the entity set.
     * \return A reference-counted pointer to the parallel communicator.
     */
    virtual Teuchos::RCP<const Teuchos::Comm<int> > communicator() const;
    //@}

    //@{
    //! Entity access functions.
    /*!
     * \brief Get an iterator over a subset of the entity set that satisfies
     * the given predicate.
     * \param entity_type The type of entity to get an iterator over.
     * \param predicate A predicate to select the entity set to iterate over.
     * \return An iterator to the entities that satisfy the predicate.
     */
    virtual AbstractIterator<Entity>
    entityIterator(
	const EntityType entity_type,
	const std::function<bool(const Entity&)>& predicate = selectAll 
	) const;

    /*!
     * \brief Given an EntityId, get the entity.
     * \param entity_id Get the entity with this id.
     * \param entity The entity with the given id.
     */
    virtual void getEntity( const EntityId entity_id, 
			    Entity& entity ) const;
    //@}

    //@{
    //! Geometric data functions.
    /*!
     * \brief Return the largest physical dimension of the entities in the
     * set.  
     * \return The physical dimension of the set.
     */
    virtual int physicalDimension() const;

    /*!
     * \brief Get the local bounding box of entities of the set.
     * \return A Cartesian box the bounds all local entities in the set.
     */
    virtual void localBoundingBox( Teuchos::Tuple<double,6>& bounds ) const;

    /*!
     * \brief Get the global bounding box of entities of the set.
     * \return A Cartesian box the bounds all global entities in the set.
     */
    virtual void globalBoundingBox( Teuchos::Tuple<double,6>& bounds ) const;
    //@}

    //@{
    //! EntitySet modification functions.
    /*!
     * \brief Indicate that the entity set will be modified.
     */
    virtual void startModification();

    /*!
     * \brief Add an entity to the set. This function can only be called if
     * entity set has been notified of modification.
     * \param entity Add this entity to the set.
     */
    virtual void addEntity( const Entity& entity );

    /*!
     * \brief Indicate that modification of the entity set is complete.
     */
    virtual void endModification();
    //@}

    /*!
     * \brief Select all entities predicate.
     */
    static inline bool selectAll( const Entity& entity )
    { return true; }
};

//---------------------------------------------------------------------------//
// AbstractBuildableObjectPolicy implementation.
//---------------------------------------------------------------------------//
template<>
class AbstractBuildableObjectPolicy<EntitySet>
{
  public:

    typedef EntitySet object_type;

    static std::string 
    objectType( const EntitySet& entity_set )
    {
	return entity_set.name();
    }

    static Teuchos::RCP<DataTransferKit::AbstractBuilder<EntitySet> > 
    getBuilder()
    {
	return EntitySet::getBuilder();
    }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ENTITYSET_HPP

//---------------------------------------------------------------------------//
// end DTK_EntitySet.hpp
//---------------------------------------------------------------------------//
