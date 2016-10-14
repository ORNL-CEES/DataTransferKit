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
 * \brief DTK_LibmeshManager.cpp
 * \author Stuart R. Slattery
 * \brief High-level manager for Libmesh.
 */
//---------------------------------------------------------------------------//

#include "DTK_LibmeshManager.hpp"
#include "DTK_LibmeshEntityIntegrationRule.hpp"
#include "DTK_LibmeshEntityLocalMap.hpp"
#include "DTK_LibmeshEntitySet.hpp"
#include "DTK_LibmeshNodalShapeFunction.hpp"
#include "DTK_LibmeshVariableField.hpp"

#include <DTK_BasicEntityPredicates.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
//! Default constructor.
LibmeshManager::LibmeshManager(
    const Teuchos::RCP<libMesh::MeshBase> &libmesh_mesh,
    const Teuchos::RCP<libMesh::System> &libmesh_system )
    : d_mesh( libmesh_mesh )
    , d_system( libmesh_system )
{
    SelectAllPredicate pred;
    buildFunctionSpace( pred.getFunction() );
}

//---------------------------------------------------------------------------//
//! Subdomain constructor.
LibmeshManager::LibmeshManager(
    const Teuchos::RCP<libMesh::MeshBase> &libmesh_mesh,
    const Teuchos::RCP<libMesh::System> &libmesh_system,
    const Teuchos::Array<libMesh::subdomain_id_type> &subdomain_ids )
    : d_mesh( libmesh_mesh )
    , d_system( libmesh_system )
{
    Teuchos::Array<int> int_ids( subdomain_ids.begin(), subdomain_ids.end() );
    BlockPredicate pred( int_ids );
    buildFunctionSpace( pred.getFunction() );
}

//---------------------------------------------------------------------------//
//! Boundary constructor.
LibmeshManager::LibmeshManager(
    const Teuchos::RCP<libMesh::MeshBase> &libmesh_mesh,
    const Teuchos::RCP<libMesh::System> &libmesh_system,
    const Teuchos::Array<libMesh::boundary_id_type> &boundary_ids )
    : d_mesh( libmesh_mesh )
    , d_system( libmesh_system )
{
    Teuchos::Array<int> int_ids( boundary_ids.begin(), boundary_ids.end() );
    BoundaryPredicate pred( int_ids );
    buildFunctionSpace( pred.getFunction() );
}

//---------------------------------------------------------------------------//
// Get the function space over which the mesh and its fields are defined.
Teuchos::RCP<FunctionSpace> LibmeshManager::functionSpace() const
{
    return d_function_space;
}

//---------------------------------------------------------------------------//
// Build the function space.
void LibmeshManager::buildFunctionSpace( const PredicateFunction &pred )
{
    Teuchos::RCP<EntitySet> entity_set =
        Teuchos::rcp( new LibmeshEntitySet( d_mesh ) );

    Teuchos::RCP<EntityLocalMap> local_map =
        Teuchos::rcp( new LibmeshEntityLocalMap( d_mesh, d_system ) );

    Teuchos::RCP<EntityShapeFunction> shape_function =
        Teuchos::rcp( new LibmeshNodalShapeFunction( d_mesh, d_system ) );

    Teuchos::RCP<EntityIntegrationRule> integration_rule =
        Teuchos::rcp( new LibmeshEntityIntegrationRule() );

    d_function_space = Teuchos::rcp( new FunctionSpace(
        entity_set, local_map, shape_function, integration_rule, pred ) );
}

//---------------------------------------------------------------------------//
// Build a field vector from a variable.
Teuchos::RCP<FieldMultiVector>
LibmeshManager::createFieldMultiVector( const std::string &variable_name )
{
    Teuchos::RCP<Field> field = Teuchos::rcp(
        new LibmeshVariableField( d_mesh, d_system, variable_name ) );
    return Teuchos::rcp(
        new FieldMultiVector( field, d_function_space->entitySet() ) );
}

//---------------------------------------------------------------------------//
// Get the entity set over which the fields are defined.
Teuchos::RCP<EntitySet> LibmeshManager::entitySet() const
{
    return d_function_space->entitySet();
}

//---------------------------------------------------------------------------//
// Get the local map for entities supporting the function.
Teuchos::RCP<EntityLocalMap> LibmeshManager::localMap() const
{
    return d_function_space->localMap();
}

//---------------------------------------------------------------------------//
// Get the shape function for entities supporting the function.
Teuchos::RCP<EntityShapeFunction> LibmeshManager::shapeFunction() const
{
    return d_function_space->shapeFunction();
}

//---------------------------------------------------------------------------//
// Get the integration rule for entities supporting the function.
Teuchos::RCP<EntityIntegrationRule> LibmeshManager::integrationRule() const
{
    return d_function_space->integrationRule();
}

//---------------------------------------------------------------------------//
// Get the selector function.
PredicateFunction LibmeshManager::selectFunction() const
{
    return d_function_space->selectFunction();
}

//---------------------------------------------------------------------------//
// Get the field for the given string key.
Teuchos::RCP<Field> LibmeshManager::field( const std::string &field_name ) const
{
    return Teuchos::rcp(
        new LibmeshVariableField( d_mesh, d_system, field_name ) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_LibmeshManager.cpp
//---------------------------------------------------------------------------//
