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
 * \brief DTK_ClassicMeshEntitySet_impl.hpp
 * \author Stuart R. Slattery
 * \brief ClassicMesh entity set implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CLASSICMESHENTITYSET_IMPL_HPP
#define DTK_CLASSICMESHENTITYSET_IMPL_HPP

#include "DTK_DBC.hpp"
#include "DTK_ClassicMeshElement.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_OrdinalTraits.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// ClassicMeshEntitySetIterator implementation.
//---------------------------------------------------------------------------//
// Default constructor.
template<class Mesh>
ClassicMeshEntitySetIterator<Mesh>::ClassicMeshEntitySetIterator()
{
    this->b_iterator_impl = NULL;
}

//---------------------------------------------------------------------------//
// Constructor.
template<class Mesh>
ClassicMeshEntitySetIterator<Mesh>::ClassicMeshEntitySetIterator(
    const Teuchos::RCP<ClassicMesh<Mesh> >& mesh,
    const PredicateFunction& predicate )
    : d_mesh( mesh )
    , d_current_block( 0 )
{
    DTK_REQUIRE( d_mesh->getNumBlocks() > 0 );
    d_element_it = d_mesh->d_element_gids.front().begin();
    if ( d_mesh->d_element_gids.front().end() == d_element_it )
    {
	moveToNextBlock();
    }
    this->b_iterator_impl = NULL;
    this->b_predicate = predicate;
}

//---------------------------------------------------------------------------//
// Copy constructor.
template<class Mesh>
ClassicMeshEntitySetIterator<Mesh>::ClassicMeshEntitySetIterator( 
    const ClassicMeshEntitySetIterator<Mesh>& rhs )
    : d_mesh( rhs.d_mesh )
    , d_current_block( rhs.d_current_block )
    , d_element_it( rhs.d_element_it )
{
    this->b_iterator_impl = NULL;
    this->b_predicate = rhs.b_predicate;
}

//---------------------------------------------------------------------------//
// Assignment operator.
template<class Mesh>
ClassicMeshEntitySetIterator<Mesh>& ClassicMeshEntitySetIterator<Mesh>::operator=( 
    const ClassicMeshEntitySetIterator<Mesh>& rhs )
{
    this->b_iterator_impl = NULL;
    this->b_predicate = rhs.b_predicate;
    if ( &rhs == this )
    {
	return *this;
    }
    d_mesh = rhs.d_mesh;
    d_current_block = rhs.d_current_block;
    d_element_it = rhs.d_element_it;
    return *this;
}

//---------------------------------------------------------------------------//
// Destructor.
template<class Mesh>
ClassicMeshEntitySetIterator<Mesh>::~ClassicMeshEntitySetIterator()
{
    this->b_iterator_impl = NULL;
}

//---------------------------------------------------------------------------//
// Pre-increment operator.
template<class Mesh>
EntityIterator& ClassicMeshEntitySetIterator<Mesh>::operator++()
{
    DTK_REQUIRE( d_mesh->getNumBlocks() > 0 );
    DTK_REQUIRE( d_current_block < d_mesh->getNumBlocks() );

    if ( d_element_it == d_mesh->d_element_gids[d_current_block].end() )
    {
	moveToNextBlock();
    }
    else
    {
	++d_element_it;
	if ( d_element_it == d_mesh->d_element_gids[d_current_block].end() )
	{
	    moveToNextBlock();
	}
    }
    
    return *this;
}

//---------------------------------------------------------------------------//
// Dereference operator.
template<class Mesh>
Entity& ClassicMeshEntitySetIterator<Mesh>::operator*(void)
{
    this->operator->();
    return d_entity;
}

//---------------------------------------------------------------------------//
// Dereference operator.
template<class Mesh>
Entity* ClassicMeshEntitySetIterator<Mesh>::operator->(void)
{
    d_entity =
	ClassicMeshElement<Mesh>( d_mesh.ptr(), *d_element_it, d_current_block );
    return &d_entity;
}

//---------------------------------------------------------------------------//
// Equal comparison operator.
template<class Mesh>
bool ClassicMeshEntitySetIterator<Mesh>::operator==( 
    const EntityIterator& rhs ) const
{ 
    const ClassicMeshEntitySetIterator<Mesh>* rhs_vec = 
	static_cast<const ClassicMeshEntitySetIterator<Mesh>*>(&rhs);
    const ClassicMeshEntitySetIterator<Mesh>* rhs_vec_impl = 
	static_cast<const ClassicMeshEntitySetIterator<Mesh>*>(
	    rhs_vec->b_iterator_impl);
    return ( rhs_vec_impl->d_element_it == d_element_it );
}

//---------------------------------------------------------------------------//
// Not equal comparison operator.
template<class Mesh>
bool ClassicMeshEntitySetIterator<Mesh>::operator!=( 
    const EntityIterator& rhs ) const
{
    const ClassicMeshEntitySetIterator<Mesh>* rhs_vec = 
	static_cast<const ClassicMeshEntitySetIterator<Mesh>*>(&rhs);
    const ClassicMeshEntitySetIterator<Mesh>* rhs_vec_impl = 
	static_cast<const ClassicMeshEntitySetIterator<Mesh>*>(
	    rhs_vec->b_iterator_impl);
    return ( rhs_vec_impl->d_element_it != d_element_it );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the beginning.
template<class Mesh>
EntityIterator ClassicMeshEntitySetIterator<Mesh>::begin() const
{ 
    return ClassicMeshEntitySetIterator<Mesh>( d_mesh , this->b_predicate );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the end.
template<class Mesh>
EntityIterator ClassicMeshEntitySetIterator<Mesh>::end() const
{
    ClassicMeshEntitySetIterator<Mesh> end_it( d_mesh, this->b_predicate );
    end_it.d_element_it = end_it.d_mesh->d_element_gids.back().end();
    return end_it;
}

//---------------------------------------------------------------------------//
// Create a clone of the iterator. We need this for the copy constructor
// and assignment operator to pass along the underlying implementation.
template<class Mesh>
EntityIterator* ClassicMeshEntitySetIterator<Mesh>::clone() const
{
    return new ClassicMeshEntitySetIterator<Mesh>(*this);
}

//---------------------------------------------------------------------------//
// Move the iterator to the valid block or to the end.
template<class Mesh>
void ClassicMeshEntitySetIterator<Mesh>::moveToNextBlock()
{
    DTK_REQUIRE( d_element_it == d_mesh->d_element_gids[d_current_block].end() );
    ++d_current_block;
    while ( d_current_block < d_mesh->getNumBlocks() )
    {
	d_element_it = d_mesh->d_element_gids[d_current_block].begin();
	if ( d_mesh->d_element_gids[d_current_block].end() == d_element_it )
	{
	    ++d_current_block;
	}
	else
	{
	    break;
	}
    }
}

//---------------------------------------------------------------------------//
// ClassicMeshEntitySet implementation.
//---------------------------------------------------------------------------//
// Constructor.
template<class Mesh>
ClassicMeshEntitySet<Mesh>::ClassicMeshEntitySet(
    const Teuchos::RCP<ClassicMesh<Mesh> >& mesh )
    : d_mesh( mesh )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the parallel communicator for the entity set.
template<class Mesh>
Teuchos::RCP<const Teuchos::Comm<int> >
ClassicMeshEntitySet<Mesh>::communicator() const
{
    return d_mesh->comm();
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entities in the set.
template<class Mesh>
int ClassicMeshEntitySet<Mesh>::physicalDimension() const
{
    return d_mesh->dim();
}

//---------------------------------------------------------------------------//
// Given an EntityId, get the entity.
template<class Mesh>
void ClassicMeshEntitySet<Mesh>::getEntity( const EntityId entity_id,
					    const int topological_dimension,
					    Entity& entity ) const
{
    entity = ClassicMeshElement<Mesh>( d_mesh.ptr(),
				       entity_id,
				       d_mesh->elementBlockId(entity_id) );
}

//---------------------------------------------------------------------------//
// Get an iterator over a subset of the entity set that satisfies the given
// predicate.
template<class Mesh>
EntityIterator ClassicMeshEntitySet<Mesh>::entityIterator(
    const int topological_dimension,
    const PredicateFunction& predicate ) const
{
    DTK_REQUIRE( d_mesh->dim() == topological_dimension );
    return ClassicMeshEntitySetIterator<Mesh>( d_mesh, predicate );
}

//---------------------------------------------------------------------------//
// Given an entity, get the entities of the given type that are adjacent to
// it.
template<class Mesh>
void ClassicMeshEntitySet<Mesh>::getAdjacentEntities(
    const Entity& entity,
    const int adjacent_dimension,
    Teuchos::Array<Entity>& adjacent_entities ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_CLASSICMESHENTITYSET_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_ClassicMeshEntitySet_impl.hpp
//---------------------------------------------------------------------------//
