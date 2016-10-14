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
 * \brief DTK_LibmeshEntityIterator_impl.hpp
 * \author Stuart R. Slattery
 * \brief Entity iterator interface.
 */
//---------------------------------------------------------------------------//

#ifndef LIBMESHDTKADAPTERS_LIBMESHENTITYITERATOR_IMPL_HPP
#define LIBMESHDTKADAPTERS_LIBMESHENTITYITERATOR_IMPL_HPP

#include <type_traits>

#include "DTK_LibmeshEntity.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Default constructor.
template <class LibmeshGeomIterator>
LibmeshEntityIterator<LibmeshGeomIterator>::LibmeshEntityIterator()
{ /* ... */
}

//---------------------------------------------------------------------------//
// Constructor.
template <class LibmeshGeomIterator>
LibmeshEntityIterator<LibmeshGeomIterator>::LibmeshEntityIterator(
    LibmeshGeomIterator libmesh_iterator,
    LibmeshGeomIterator libmesh_iterator_begin,
    LibmeshGeomIterator libmesh_iterator_end,
    const Teuchos::Ptr<libMesh::MeshBase> &libmesh_mesh,
    const Teuchos::Ptr<LibmeshAdjacencies> &adjacencies,
    const std::function<bool( Entity )> &predicate )
    : d_libmesh_iterator( libmesh_iterator )
    , d_libmesh_iterator_begin( libmesh_iterator_begin )
    , d_libmesh_iterator_end( libmesh_iterator_end )
    , d_libmesh_mesh( libmesh_mesh )
    , d_adjacencies( adjacencies )
{
    this->b_predicate = predicate;
}

//---------------------------------------------------------------------------//
// Copy constructor.
template <class LibmeshGeomIterator>
LibmeshEntityIterator<LibmeshGeomIterator>::LibmeshEntityIterator(
    const LibmeshEntityIterator<LibmeshGeomIterator> &rhs )
    : d_libmesh_iterator( rhs.d_libmesh_iterator )
    , d_libmesh_iterator_begin( rhs.d_libmesh_iterator_begin )
    , d_libmesh_iterator_end( rhs.d_libmesh_iterator_end )
    , d_libmesh_mesh( rhs.d_libmesh_mesh )
    , d_adjacencies( rhs.d_adjacencies )
{
    this->b_predicate = rhs.b_predicate;
}

//---------------------------------------------------------------------------//
// Assignment operator.
template <class LibmeshGeomIterator>
LibmeshEntityIterator<LibmeshGeomIterator> &
LibmeshEntityIterator<LibmeshGeomIterator>::
operator=( const LibmeshEntityIterator<LibmeshGeomIterator> &rhs )
{
    this->b_predicate = rhs.b_predicate;
    if ( &rhs == this )
    {
        return *this;
    }
    d_libmesh_iterator = rhs.d_libmesh_iterator;
    d_libmesh_iterator_begin = rhs.d_libmesh_iterator_begin;
    d_libmesh_iterator_end = rhs.d_libmesh_iterator_end;
    d_libmesh_mesh = rhs.d_libmesh_mesh;
    d_adjacencies = rhs.d_adjacencies;
    return *this;
}

//---------------------------------------------------------------------------//
// Pre-increment operator.
template <class LibmeshGeomIterator>
EntityIterator &LibmeshEntityIterator<LibmeshGeomIterator>::operator++()
{
    ++d_libmesh_iterator;
    return *this;
}

//---------------------------------------------------------------------------//
// Dereference operator.
template <class LibmeshGeomIterator>
Entity &LibmeshEntityIterator<LibmeshGeomIterator>::operator*( void )
{
    this->operator->();
    return d_current_entity;
}

//---------------------------------------------------------------------------//
// Dereference operator.
template <class LibmeshGeomIterator>
Entity *LibmeshEntityIterator<LibmeshGeomIterator>::operator->( void )
{
    d_current_entity = LibmeshEntity<typename std::remove_pointer<
        typename LibmeshGeomIterator::value_type>::type>(
        Teuchos::ptr( *d_libmesh_iterator ), d_libmesh_mesh, d_adjacencies );
    return &d_current_entity;
}

//---------------------------------------------------------------------------//
// Equal comparison operator.
template <class LibmeshGeomIterator>
bool LibmeshEntityIterator<LibmeshGeomIterator>::
operator==( const EntityIterator &rhs ) const
{
    const LibmeshEntityIterator *rhs_it =
        static_cast<const LibmeshEntityIterator *>( &rhs );
    const LibmeshEntityIterator *rhs_it_impl =
        static_cast<const LibmeshEntityIterator *>(
            rhs_it->b_iterator_impl.get() );
    return ( rhs_it_impl->d_libmesh_iterator == d_libmesh_iterator );
}

//---------------------------------------------------------------------------//
// Not equal comparison operator.
template <class LibmeshGeomIterator>
bool LibmeshEntityIterator<LibmeshGeomIterator>::
operator!=( const EntityIterator &rhs ) const
{
    const LibmeshEntityIterator *rhs_it =
        static_cast<const LibmeshEntityIterator *>( &rhs );
    const LibmeshEntityIterator *rhs_it_impl =
        static_cast<const LibmeshEntityIterator *>(
            rhs_it->b_iterator_impl.get() );
    return ( rhs_it_impl->d_libmesh_iterator != d_libmesh_iterator );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the beginning.
template <class LibmeshGeomIterator>
EntityIterator LibmeshEntityIterator<LibmeshGeomIterator>::begin() const
{
    return LibmeshEntityIterator( d_libmesh_iterator_begin,
                                  d_libmesh_iterator_begin,
                                  d_libmesh_iterator_end, d_libmesh_mesh,
                                  d_adjacencies, this->b_predicate );
}

//---------------------------------------------------------------------------//
// An iterator assigned to the end.
template <class LibmeshGeomIterator>
EntityIterator LibmeshEntityIterator<LibmeshGeomIterator>::end() const
{
    return LibmeshEntityIterator( d_libmesh_iterator_end,
                                  d_libmesh_iterator_begin,
                                  d_libmesh_iterator_end, d_libmesh_mesh,
                                  d_adjacencies, this->b_predicate );
}

//---------------------------------------------------------------------------//
// Create a clone of the iterator. We need this for the copy constructor
// and assignment operator to pass along the underlying implementation.
template <class LibmeshGeomIterator>
std::unique_ptr<EntityIterator>
LibmeshEntityIterator<LibmeshGeomIterator>::clone() const
{
    return std::unique_ptr<EntityIterator>(
        new LibmeshEntityIterator( *this ) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end LIBMESHDTKADAPTERS_LIBMESHENTITYITERATOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_LibmeshEntityIterator_impl.hpp
//---------------------------------------------------------------------------//
