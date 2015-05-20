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
 * \brief DTK_ClassicMeshEntitySet.hpp
 * \author Stuart R. Slattery
 * \brief ClassicMesh entity set implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CLASSICMESHENTITYSET_HPP
#define DTK_CLASSICMESHENTITYSET_HPP

#include <unordered_map>

#include "DTK_Types.hpp"
#include "DTK_EntitySet.hpp"
#include "DTK_Entity.hpp"
#include "DTK_EntityIterator.hpp"
#include "DTK_ClassicMesh.hpp"
#include "DTK_Classic_MeshTraits.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class ClassicMeshEntitySetIterator
  \brief Implementation of iterator over entities in a classicmesh set.
*/
template<class Mesh>
class ClassicMeshEntitySetIterator : public EntityIterator
{
  public:

    // Typedefs.
    typedef Classic::MeshTraits<Mesh> MT;
    typedef typename MT::global_ordinal_type GlobalOrdinal;

    // Default constructor.
    ClassicMeshEntitySetIterator();

    // Constructor.
    ClassicMeshEntitySetIterator( const Teuchos::RCP<ClassicMesh<Mesh> >& mesh,
				  const PredicateFunction& predicate );

    // Copy constructor.
    ClassicMeshEntitySetIterator( const ClassicMeshEntitySetIterator<Mesh>& rhs );

    // Destructor.
    ~ClassicMeshEntitySetIterator();
    
    /*!
     * \brief Assignment operator.
     */
    ClassicMeshEntitySetIterator<Mesh>&
    operator=( const ClassicMeshEntitySetIterator<Mesh>& rhs );

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

    // Move the iterator to the valid block or to the end.
    void moveToNextBlock();
    
  private:

    // Classic mesh.
    Teuchos::RCP<ClassicMesh<Mesh> > d_mesh;

    // Current block id.
    int d_current_block;

    // Iterator over the current block.
    typename Teuchos::ArrayRCP<const GlobalOrdinal>::const_iterator d_element_it;

    // Pointer to the current entity.
    Entity d_entity;
};

//---------------------------------------------------------------------------//
/*!
  \class ClassicMeshEntitySet
  \brief ClassicMesh implementation of the entity set interface.
*/
//---------------------------------------------------------------------------//
template<class Mesh>
class ClassicMeshEntitySet : public EntitySet
{
  public:

    /*!
     * \brief Constructor.
     */
    ClassicMeshEntitySet( const Teuchos::RCP<ClassicMesh<Mesh> >& mesh );

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
     * \param topological_dimension Get the entity with this topological
     * dimension.
     * \param entity The entity with the given id.
     */
    void getEntity( const EntityId entity_id,
		    const int topological_dimension,
		    Entity& entity ) const override;

    /*!
     * \brief Get an iterator over a subset of the entity set that satisfies
     * the given predicate.
     * \param topological_dimension The topological dimension of entity to get
     * an iterator for.
     * \param predicate A predicate to select the entity set to iterate over.
     * \return An iterator to the entities that satisfy the predicate.
     */
    EntityIterator entityIterator(
	const int topological_dimension,
	const PredicateFunction& predicate ) const override;

    /*!
     * \brief Given an entity, get the entities of the given type that are
     * adjacent to it.
     */
    virtual void getAdjacentEntities(
	const Entity& entity,
	const int adjacent_dimension,
	Teuchos::Array<Entity>& adjacent_entities ) const override;
    //@}

  private:

    // Classic mesh.
    Teuchos::RCP<ClassicMesh<Mesh> > d_mesh;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_ClassicMeshEntitySet_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_CLASSICMESHENTITYSET_HPP

//---------------------------------------------------------------------------//
// end DTK_ClassicMeshEntitySet.hpp
//---------------------------------------------------------------------------//
