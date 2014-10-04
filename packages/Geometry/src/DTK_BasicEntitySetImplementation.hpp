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
 * \brief DTK_BasicEntitySetImplementation.hpp
 * \author Stuart R. Slattery
 * \brief Basic entity set implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BASICENTITYSETIMPLEMENTATION_HPP
#define DTK_BASICENTITYSETIMPLEMENTATION_HPP

#include <unordered_map>

#include "DTK_EntitySet.hpp"
#include "DTK_GeometricEntity.hpp"
#include "DTK_Box.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class BasicEntitySetImplementation
  \brief Basic implementation of the entity set interface.
*/
//---------------------------------------------------------------------------//
class BasicEntitySetImplementation : public EntitySet
{
  public:

    //! Entity map typedefs.
    typedef std::unordered_map<EntityId,Teuchos::RCP<GeometricEntity> >::iterator
    EntityIterator;
    typedef std::unordered_map<EntityId,Teuchos::RCP<GeometricEntity> >::const_iterator
    ConstEntityIterator;
    typedef std::pair<EntityId,Teuchos::RCP<GeometricEntity> > EntityIdPair;

    /*!
     * \brief Default constructor.
     */
    BasicEntitySetImplementation();

    /*!
     * \brief Constructor.
     */
    BasicEntitySetImplementation( 
	const Teuchos::RCP<const Teuchos::Comm<int> > comm,
	const int physical_dimension );

    /*!
     * \brief Destructor.
     */
    ~BasicEntitySetImplementation();

    //@{
    //! Identification functions.
    /*!
     * \brief Assign a parallel communicator to the entity set. This will only
     * be done immediately after construct through the AbstractBuilder
     * interface.
     * \param comm The communicator to set with the entity set.
     */
    void assignCommunicator( 
	const Teuchos::RCP<const Teuchos::Comm<int> >& comm );
    /*!
     * \brief Return a string indicating the derived entity set type.
     * \return A string key for the derived set type.
     */
    std::string entitySetType() const;
    //@}

    //@{
    //! Parallel functions.
    /*!
     * \brief Get the parallel communicator for the entity set.
     * \return A reference-counted pointer to the parallel communicator.
     */
    Teuchos::RCP<const Teuchos::Comm<int> > communicator() const;
    //@}

    //@{
    //! Entity access functions.
    /*!
     * \brief Get the local number of entities in the set of the given
     * parametric dimension.
     * \param parametric_dimension Get the number of entities of this
     * parametric dimension.
     * \return The local number of entities in the set.
     */
    std::size_t localNumberOfEntities( 
	const int parametric_dimension ) const;

    /*!
     * \brief Get the global number of entities in the set of the given
     * parametric dimension.
     * \param parametric_dimension Get the number of entities of this
     * parametric dimension.
     * \return The global number of entities in the set.
     */
    std::size_t globalNumberOfEntities(
	const int parametric_dimension ) const;
    
    /*!
     * \brief Get the identifiers for all local entities in the set of a given
     * parametric dimension.
     * \param parametric_dimension Get the entities of this parametric
     * dimension.
     * \param entity_ids A view into an array of size localNumberOfEntities(
     * parametric_dimension ). Write the entity ids into this array.
     */
    void localEntityIds( 
	const int parametric_dimension,
	const Teuchos::ArrayView<EntityId>& entity_ids ) const;

    /*!
     * \brief Given an EntityId, get the entity.
     * \param entity_id Get the entity with this id.
     * \param entity The entity with the given id.
     */
    void getEntity( const EntityId entity_id, 
		    Teuchos::RCP<GeometricEntity>& entity ) const;
    //@}

    //@{
    //! Geometric data functions.
    /*!
     * \brief Return the physical dimension of the entities in the set.
     * \return The physical dimension of the set.
     */
    int physicalDimension() const;

    /*!
     * \brief Get the local bounding box of entities of the set.
     * \return A Cartesian box the bounds all local entities in the set.
     */
    void localBoundingBox( Box& bounding_box ) const;

    /*!
     * \brief Get the global bounding box of entities of the set.
     * \return A Cartesian box the bounds all global entities in the set.
     */
    void globalBoundingBox( Box& bounding_box ) const;
    //@}

    //@{
    //! BasicEntitySetImplementation modification functions.
    /*!
     * \brief Indicate that the entity set will be modified.
     */
    void startModification();

    /*!
     * \brief Add an entity to the set. This function can only be called if
     * entity set has been notified of modification.
     * \param entity Add this entity to the set.
     */
    void addEntity( const Teuchos::RCP<GeometricEntity>& entity );

    /*!
     * \brief Indicate that modification of the entity set is complete.
     */
    void endModification();
    //@}

  private:

    // Parallel communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > d_comm;

    // Physical dimension.
    int d_dimension;

    // Entity-to-id map.
    std::unordered_map<EntityId,Teuchos::RCP<GeometricEntity> > d_entities;
 };

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_BASICENTITYSETIMPLEMENTATION_HPP

//---------------------------------------------------------------------------//
// end DTK_BasicEntitySetImplementation.hpp
//---------------------------------------------------------------------------//
