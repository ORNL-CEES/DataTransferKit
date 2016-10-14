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
 * \brief DTK_LibmeshManager.hpp
 * \author Stuart R. Slattery
 * \brief High-level manager for Libmesh.
 */
//---------------------------------------------------------------------------//

#ifndef LIBMESHDTKADAPTERS_LIBMESHMANAGER_HPP
#define LIBMESHDTKADAPTERS_LIBMESHMANAGER_HPP

#include <functional>
#include <string>

#include "DTK_ClientManager.hpp"
#include "DTK_FieldMultiVector.hpp"
#include "DTK_FunctionSpace.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_RCP.hpp>

#include <libmesh/mesh_base.h>
#include <libmesh/system.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class LibmeshManager
  \brief High-level manager for Libmesh mesh.

  This manager provides a high-level class for automated construction of DTK
  interface objects. A user is not required to use this class but rather could
  use it to reduce code for certain implementations.
*/
//---------------------------------------------------------------------------//
class LibmeshManager : public ClientManager
{
  public:
    /*!
     * \brief Default constructor.
     *
     * \param libmesh_mesh Libmesh mesh.
     *
     * \param libsystem_system Libsystem system.
     *
     * \param entity_type The type of entities in the mesh that will be
     * mapped.
     */
    LibmeshManager( const Teuchos::RCP<libMesh::MeshBase> &libmesh_mesh,
                    const Teuchos::RCP<libMesh::System> &libmesh_system );

    /*!
     * \brief Subdomain constructor.
     *
     * \param libmesh_mesh Libmesh mesh.
     *
     * \param libsystem_system Libsystem system.
     *
     * \param subdomain_ids The subdomain ids to map.
     *
     * \param entity_type The type of entities in the mesh that will be
     * mapped.
     */
    LibmeshManager(
        const Teuchos::RCP<libMesh::MeshBase> &libmesh_mesh,
        const Teuchos::RCP<libMesh::System> &libmesh_system,
        const Teuchos::Array<libMesh::subdomain_id_type> &subdomain_ids );

    /*!
     * \brief Boundary constructor.
     *
     * \param libmesh_mesh Libmesh mesh.
     *
     * \param libsystem_system Libsystem system.
     *
     * \param boundary_ids The boundary ids to map.
     *
     * \param entity_type The type of entities in the mesh that will be
     * mapped.
     */
    LibmeshManager(
        const Teuchos::RCP<libMesh::MeshBase> &libmesh_mesh,
        const Teuchos::RCP<libMesh::System> &libmesh_system,
        const Teuchos::Array<libMesh::boundary_id_type> &boundary_ids );

    /*!
     * \brief Get the function space over which the mesh and its fields are
     * defined.
     */
    Teuchos::RCP<FunctionSpace> functionSpace() const;

    /*!
     * \brief Given a variable name, build a field vector.
     */
    Teuchos::RCP<FieldMultiVector>
    createFieldMultiVector( const std::string &variable_name );

    //@{
    //! ClientManager interface implementation.
    /*!
     * \brief Get the entity set over which the fields are defined.
     */
    Teuchos::RCP<EntitySet> entitySet() const override;

    /*!
     * \brief Get the local map for entities supporting the function.
     */
    Teuchos::RCP<EntityLocalMap> localMap() const override;

    /*!
     * \brief Get the shape function for entities supporting the function.
     */
    Teuchos::RCP<EntityShapeFunction> shapeFunction() const override;

    /*!
     * \brief Get the integration rule for entities supporting the function.
     */
    Teuchos::RCP<EntityIntegrationRule> integrationRule() const override;

    /*!
     * \brief Get the selector function.
     */
    PredicateFunction selectFunction() const override;

    /*!
     * \brief Get the field for the given string key.
     */
    Teuchos::RCP<Field> field( const std::string &field_name ) const override;
    //@}

  private:
    // Build the function space.
    void buildFunctionSpace( const PredicateFunction &pred );

  private:
    // The mesh.
    Teuchos::RCP<libMesh::MeshBase> d_mesh;

    // The system.
    Teuchos::RCP<libMesh::System> d_system;

    // The function space over which the mesh and its fields are defined.
    Teuchos::RCP<FunctionSpace> d_function_space;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end LIBMESHDTKADAPTERS_LIBMESHMANAGER_HPP

//---------------------------------------------------------------------------//
// end DTK_LibmeshManager.hpp
//---------------------------------------------------------------------------//
