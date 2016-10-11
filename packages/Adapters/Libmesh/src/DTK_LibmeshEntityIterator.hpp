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
 * \brief DTK_LibmeshEntityIterator.hpp
 * \author Stuart R. Slattery
 * \brief Entity iterator interface.
 */
//---------------------------------------------------------------------------//

#ifndef LIBMESHDTKADAPTERS_LIBMESHENTITYITERATOR_HPP
#define LIBMESHDTKADAPTERS_LIBMESHENTITYITERATOR_HPP

#include <functional>
#include <vector>

#include "DTK_LibmeshAdjacencies.hpp"

#include <DTK_Entity.hpp>
#include <DTK_EntityIterator.hpp>

#include <Teuchos_Ptr.hpp>

#include <libmesh/elem.h>
#include <libmesh/mesh_base.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class LibmeshEntityIterator
  \brief Libmesh mesh entity iterator implementation.
*/
//---------------------------------------------------------------------------//
template <class LibmeshGeomIterator>
class LibmeshEntityIterator : public EntityIterator
{
  public:
    /*!
     * \brief Default constructor.
     */
    LibmeshEntityIterator();

    /*!
     * \brief Constructor.
     */
    LibmeshEntityIterator( LibmeshGeomIterator libmesh_iterator,
                           LibmeshGeomIterator libmesh_iterator_begin,
                           LibmeshGeomIterator libmesh_iterator_end,
                           const Teuchos::Ptr<libMesh::MeshBase> &libmesh_mesh,
                           const Teuchos::Ptr<LibmeshAdjacencies> &adjacencies,
                           const PredicateFunction &predicate );

    /*!
     * \brief Copy constructor.
     */
    LibmeshEntityIterator(
        const LibmeshEntityIterator<LibmeshGeomIterator> &rhs );

    /*!
     * \brief Assignment operator.
     */
    LibmeshEntityIterator<LibmeshGeomIterator> &
    operator=( const LibmeshEntityIterator<LibmeshGeomIterator> &rhs );

    // Pre-increment operator.
    EntityIterator &operator++() override;

    // Dereference operator.
    Entity &operator*(void)override;

    // Dereference operator.
    Entity *operator->(void)override;

    // Equal comparison operator.
    bool operator==( const EntityIterator &rhs ) const override;

    // Not equal comparison operator.
    bool operator!=( const EntityIterator &rhs ) const override;

    // An iterator assigned to the first valid element in the iterator.
    EntityIterator begin() const override;

    // An iterator assigned to the end of all elements under the iterator.
    EntityIterator end() const override;

    // Create a clone of the iterator. We need this for the copy constructor
    // and assignment operator to pass along the underlying implementation.
    std::unique_ptr<EntityIterator> clone() const override;

  private:
    // Libmesh iterator.
    LibmeshGeomIterator d_libmesh_iterator;

    // Libmesh iterator to the beginning of the range.
    LibmeshGeomIterator d_libmesh_iterator_begin;

    // Libmesh iterator to the end of the range.
    LibmeshGeomIterator d_libmesh_iterator_end;

    // The mesh owning the entities.
    Teuchos::Ptr<libMesh::MeshBase> d_libmesh_mesh;

    // Mesh adjacencies.
    Teuchos::Ptr<LibmeshAdjacencies> d_adjacencies;

    // Current entity.
    Entity d_current_entity;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_LibmeshEntityIterator_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end LIBMESHDTKADAPTERS_LIBMESHENTITYITERATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_LibmeshEntityIterator.hpp
//---------------------------------------------------------------------------//
