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
 * \brief DTK_BasicEntitySet.hpp
 * \author Stuart R. Slattery
 * \brief Basic entity set implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BASICENTITYSET_HPP
#define DTK_BASICENTITYSET_HPP

#include <unordered_map>

#include "DTK_EntitySet.hpp"
#include "DTK_GeometricEntity.hpp"
#include "DTK_DerivedObjectRegistry.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class BasicEntitySetIterator
  \brief Iterator over entities in a basic set.
*/
class BasicEntitySetIterator : public EntitySetIterator
{
  public:

    // Default constructor.
    BasicEntitySetIterator() { /* ... */ }

    // Constructor.
    BasicEntitySetIterator(
	const std::unordered_map<EntityId,GeometricEntity>::const_iterator& 
	map_it )
	: d_map_it( map_it )
    { /* ... */ }

    // Pre-increment operator.
    EntitySetIterator& operator++()
    { 
	++d_map_it; 
	return *this; 
    }

    // Post-increment operator.
    EntitySetIterator operator++(int)
    { 
	BasicEntitySetIterator inc_it(d_map_it);
	++d_map_it;
	return inc_it;
    }

    // Dereference operator.
    const GeometricEntity& operator*(void)
    { return d_map_it->second; }

    // Dereference operator.
    const GeometricEntity* operator->(void)
    { return &(d_map_it->second); }

    // Equal comparison operator.
    bool operator==( const EntitySetIterator& rhs ) const
    { return d_map_it == 
	    Teuchos::rcp_dynamic_cast<BasicEntitySetIterator>( 
		rhs.b_iterator_impl)->d_map_it; }

    // Not equal comparison operator.
    bool operator!=( const EntitySetIterator& rhs ) const
    { return d_map_it != 
	    Teuchos::rcp_dynamic_cast<BasicEntitySetIterator>( 
		rhs.b_iterator_impl)->d_map_it; }

  private:

    // Iterator over the entity map.
    std::unordered_map<EntityId,GeometricEntity>::const_iterator d_map_it;
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
    typedef std::pair<EntityId,GeometricEntity> EntityIdPair;
    //@}

    /*!
     * \brief Default constructor.
     */
    BasicEntitySet();

    /*!
     * \brief Constructor.
     */
    BasicEntitySet( 
	const Teuchos::RCP<const Teuchos::Comm<int> > comm,
	const int physical_dimension );

    /*!
     * \brief Destructor.
     */
    ~BasicEntitySet();

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
     * \brief Get a forward iterator assigned to the beginning of the
     * entities in the set of the given parametric dimension.  
     * \param parametric_dimension Get a forward iterator to the beginning of
     * the entities of this dimension.  
     * \return a forward iterator assigned to the beginning of the entities in
     * the set of the given parametric dimension. 
     */
    std::iterator<std::forward_iterator_tag,GeometricEntity>
    entityIteratorBegin( const int parametric_dimension ) const;

    /*!
     * \brief Get a forward iterator assigned to the end of the entities in
     * the set of the given parametric dimension.
     * \param parametric_dimension Get a forward iterator to the end of the
     * entities of this dimension.
     * \return a forward iterator assigned to the end of the entities in the
     * set of the given parametric dimension.
     */
    std::iterator<std::forward_iterator_tag,GeometricEntity>
    entityIteratorEnd( const int parametric_dimension ) const;

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
		    GeometricEntity& entity ) const;
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
    void localBoundingBox( Teuchos::Tuple<double,6>& bounds ) const;

    /*!
     * \brief Get the global bounding box of entities of the set.
     * \return A Cartesian box the bounds all global entities in the set.
     */
    void globalBoundingBox( Teuchos::Tuple<double,6>& bounds ) const;
    //@}

    //@{
    //! BasicEntitySet modification functions.
    /*!
     * \brief Indicate that the entity set will be modified.
     */
    void startModification();

    /*!
     * \brief Add an entity to the set. This function can only be called if
     * entity set has been notified of modification.
     * \param entity Add this entity to the set.
     */
    void addEntity( const GeometricEntity& entity );

    /*!
     * \brief Indicate that modification of the entity set is complete.
     */
    void endModification();
    //@}

  private:

    // Given a parametric dimension, get an id-to-entity map.
    std::unordered_map<EntityId,GeometricEntity>& getEntityMap(
	const int parametric_dimension );
    
  private:

    // Parallel communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > d_comm;

    // Physical dimension.
    int d_dimension;

    // Id-to-dimension map.
    std::unordered_map<EntityId,int> d_entity_dims;

    // Id-to-entity maps.
    Teuchos::Array<std::unordered_map<EntityId,GeometricEntity> > d_entities;
};

//---------------------------------------------------------------------------//
// DerivedObjectRegistrationPolicy implementation.
//---------------------------------------------------------------------------//
template<>
class DerivedObjectRegistrationPolicy<BasicEntitySet>
{
  public:

    //! Base class type.
    typedef BasicEntitySet object_type;

    /*!
     * \brief Register a derived class with a base class.
     */
    static void registerDerivedClassWithBaseClass()
    {
	// Register the constructor with the base class
	// AbstractBuildableObject interface.
	EntitySet::setDerivedClassFactory<object_type>();
    }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_BASICENTITYSET_HPP

//---------------------------------------------------------------------------//
// end DTK_BasicEntitySet.hpp
//---------------------------------------------------------------------------//
