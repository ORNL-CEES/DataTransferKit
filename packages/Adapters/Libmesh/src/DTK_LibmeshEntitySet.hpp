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
 * \brief LibmeshDTKAdpaters_LibmeshEntitySet.hpp
 * \author Stuart R. Slattery
 * \brief Libmesh mesh entity set.
 */
//---------------------------------------------------------------------------//

#ifndef LIBMESHDTKADAPTERS_LIBMESHENTITYSET_HPP
#define LIBMESHDTKADAPTERS_LIBMESHENTITYSET_HPP

#include <functional>
#include <unordered_map>

#include "DTK_LibmeshAdjacencies.hpp"
#include "DTK_LibmeshEntity.hpp"
#include "DTK_LibmeshEntityExtraData.hpp"

#include <DTK_Entity.hpp>
#include <DTK_EntityIterator.hpp>
#include <DTK_EntitySet.hpp>
#include <DTK_Types.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>

#include <libmesh/mesh_base.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class LibmeshEntitySet
  \brief Libmesh entity set.

  Entity set implementation for Libmesh.
*/
//---------------------------------------------------------------------------//
class LibmeshEntitySet : public DataTransferKit::EntitySet
{
  public:
    /*!
     * \brief Constructor.
     */
    LibmeshEntitySet( const Teuchos::RCP<libMesh::MeshBase> &libmesh_mesh );

    //@{
    //! Parallel functions.
    /*!
     * \brief Get the parallel communicator for the entity set.
     * \return A reference-counted pointer to the parallel communicator.
     */
    Teuchos::RCP<const Teuchos::Comm<int>> communicator() const override;
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
     * \param entity The entity with the given id.
     */
    void getEntity( const DataTransferKit::EntityId entity_id,
                    const int topological_dimension,
                    DataTransferKit::Entity &entity ) const override;

    /*!
     * \brief Get a iterator of the given entity type that satisfy the given
     * predicate.
     * \param entity_type The type of entity to get a iterator for.
     * \param predicate The selection predicate.
     * \return A iterator of entities of the given type.
     */
    DataTransferKit::EntityIterator
    entityIterator( const int topological_dimension,
                    const std::function<bool( DataTransferKit::Entity )>
                        &predicate ) const override;

    /*!
     * \brief Given an entity, get the entities of the given type that are
     * adjacent to it.
     */
    void getAdjacentEntities( const DataTransferKit::Entity &entity,
                              const int adjacent_dimension,
                              Teuchos::Array<DataTransferKit::Entity>
                                  &adjacent_entities ) const override;

    /*!
     * \brief Provide a one line description of the object.
     */
    std::string description() const override
    {
        return std::string( "libMesh Mesh" );
    }
    //@}

  private:
    // Get adjacent entities implementation.
    template <class FromGeomType, class ToGeomType>
    void getAdjacentEntitiesImpl(
        const DataTransferKit::Entity &entity,
        Teuchos::Array<DataTransferKit::Entity> &adjacent_entities ) const;

  private:
    // Libmesh mesh.
    Teuchos::RCP<libMesh::MeshBase> d_libmesh_mesh;

    // Mesh adjacencies.
    Teuchos::RCP<LibmeshAdjacencies> d_adjacencies;
};

//---------------------------------------------------------------------------//
// Template functions.
//---------------------------------------------------------------------------//
// Get adjacent entities implementation.
template <class FromGeomType, class ToGeomType>
void LibmeshEntitySet::getAdjacentEntitiesImpl(
    const DataTransferKit::Entity &entity,
    Teuchos::Array<DataTransferKit::Entity> &adjacent_entities ) const
{
    Teuchos::Array<Teuchos::Ptr<ToGeomType>> adjacent_libmesh;
    d_adjacencies->getLibmeshAdjacencies(
        Teuchos::rcp_dynamic_cast<LibmeshEntityExtraData<FromGeomType>>(
            entity.extraData() )
            ->d_libmesh_geom,
        adjacent_libmesh );

    adjacent_entities.resize( adjacent_libmesh.size() );
    typename Teuchos::Array<Teuchos::Ptr<ToGeomType>>::iterator libmesh_it;
    Teuchos::Array<DataTransferKit::Entity>::iterator dtk_it;
    for ( libmesh_it = adjacent_libmesh.begin(),
          dtk_it = adjacent_entities.begin();
          libmesh_it != adjacent_libmesh.end(); ++libmesh_it, ++dtk_it )
    {
        *dtk_it = LibmeshEntity<ToGeomType>( *libmesh_it, d_libmesh_mesh.ptr(),
                                             d_adjacencies.ptr() );
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end LIBMESHDTKADAPTERS_LIBMESHENTITYSET_HPP

//---------------------------------------------------------------------------//
// end DTK_LibmeshEntitySet.hpp
//---------------------------------------------------------------------------//
