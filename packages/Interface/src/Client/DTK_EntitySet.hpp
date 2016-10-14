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
 * \brief DTK_EntitySet.hpp
 * \author Stuart R. Slattery
 * \brief Geometric entity set interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYSET_HPP
#define DTK_ENTITYSET_HPP

#include <functional>

#include "DTK_Entity.hpp"
#include "DTK_EntityIterator.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_RCP.hpp>
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
class EntitySet : public Teuchos::Describable
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
    //! Parallel functions.
    /*!
     * \brief Get the parallel communicator for the entity set.
     *
     * \return A reference-counted pointer to the parallel communicator.
     */
    virtual Teuchos::RCP<const Teuchos::Comm<int>> communicator() const = 0;
    //@}

    //@{
    //! Geometric data functions.
    /*!
     * \brief Return the largest physical dimension of the entities in the
     * set.
     *
     * \return The physical dimension of the set.
     */
    virtual int physicalDimension() const = 0;

    /*!
     * \brief Get the local bounding box of entities of the set.
     *
     * \return A Cartesian box the bounds all local entities in the set.
     */
    virtual void localBoundingBox( Teuchos::Tuple<double, 6> &bounds ) const;

    /*!
     * \brief Get the global bounding box of entities of the set.
     *
     * Requires global communication: a single all-reduce call
     *
     * \return A Cartesian box the bounds all global entities in the set.
     */
    virtual void globalBoundingBox( Teuchos::Tuple<double, 6> &bounds ) const;
    //@}

    //@{
    //! Entity access functions.
    /*!
     * \brief Given an EntityId and topological dimension, get the entity.
     *
     * \param entity_id Get the entity with this id.
     *
     * \param topological_dimension Get the entity with this topological
     * dimension.
     *
     * \param entity The entity with the given id.
     */
    virtual void getEntity( const EntityId entity_id,
                            const int topological_dimension,
                            Entity &entity ) const = 0;

    /*!
     * \brief Get an iterator of the given entity type that satisfies the
     * given predicate.
     *
     * \param topological_dimension The topological dimension of entity to get
     * an iterator for.
     *
     * \param predicate The selection predicate.
     *
     * \return A iterator of entities of the given type.
     */
    virtual EntityIterator
    entityIterator( const int topological_dimension,
                    const PredicateFunction &predicate = selectAll ) const = 0;

    /*!
     * \brief Given an entity, get the entities of the given type that are
     * adjacent to it.
     *
     * \param entity Get the adjacencies of this entity.
     *
     * \param adjacent_dimension Get adjacencies of this topological
     * dimension.
     */
    virtual void
    getAdjacentEntities( const Entity &entity, const int adjacent_dimension,
                         Teuchos::Array<Entity> &adjacent_entities ) const = 0;
    //@}

    //@{
    //! Teuchos::Describable interface.
    /*!
     * \brief Provide a one line description of the object.
     */
    virtual std::string description() const override;

    /*!
     * \brief Provide a verbose description of the object.
     */
    virtual void
    describe( Teuchos::FancyOStream &out,
              const Teuchos::EVerbosityLevel verb_level =
                  Teuchos::Describable::verbLevel_default ) const override;
    //@}

    /*!
     * \brief Select all entities predicate.
     */
    static inline bool selectAll( Entity ) { return true; }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ENTITYSET_HPP

//---------------------------------------------------------------------------//
// end DTK_EntitySet.hpp
//---------------------------------------------------------------------------//
