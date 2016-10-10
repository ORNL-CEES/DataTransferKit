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
 * \brief DTK_MoabEntitySet.hpp
 * \author Stuart R. Slattery
 * \brief Moab mesh entity set.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOABENTITYSET_HPP
#define DTK_MOABENTITYSET_HPP

#include <functional>
#include <unordered_set>

#include "DTK_EntitySet.hpp"
#include "DTK_Types.hpp"
#include "DTK_Entity.hpp"
#include "DTK_EntityIterator.hpp"
#include "DTK_MoabMeshSetIndexer.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>

#include <moab/ParallelComm.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class MoabEntitySet
  \brief Moab entity set.

  Entity set implementation for Moab.
*/
//---------------------------------------------------------------------------//
class MoabEntitySet : public EntitySet
{
  public:

    /*!
     * \brief Constructor.
     */
    MoabEntitySet( const Teuchos::RCP<moab::ParallelComm>& moab_mesh,
                   const Teuchos::RCP<MoabMeshSetIndexer> set_indexer );

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
     * \brief Return the largest physical dimension of the entities in the
     * set.
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
     * \brief Get a iterator of the given entity type that satisfy the given
     * predicate.
     * \param topological_dimension The topological dimension of entity to get
     * an iterator for.
     * \param predicate The selection predicate.
     * \return A iterator of entities of the given type.
     */
    EntityIterator entityIterator(
        const int topological_dimension,
        const PredicateFunction& predicate ) const override;

    /*!
     * \brief Given an entity, get the entities of the given type that are
     * adjacent to it.
     */
    void getAdjacentEntities(
        const Entity& entity,
        const int adjacent_dimension,
        Teuchos::Array<Entity>& adjacent_entities ) const override;

    /*!
     * \brief Provide a one line description of the object.
     */
    std::string description() const
    { return std::string("Moab Mesh"); }
    //@}

  private:

    // Moab mesh.
    Teuchos::RCP<moab::ParallelComm> d_moab_mesh;

    // Mesh set indexer.
    Teuchos::RCP<MoabMeshSetIndexer> d_set_indexer;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MOABENTITYSET_HPP

//---------------------------------------------------------------------------//
// end DTK_MoabEntitySet.hpp
//---------------------------------------------------------------------------//
