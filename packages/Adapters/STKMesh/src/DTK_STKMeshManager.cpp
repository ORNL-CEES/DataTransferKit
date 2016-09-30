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
 * \brief DTK_STKMeshManager.cpp
 * \author Stuart R. Slattery
 * \brief High-level manager for STK mesh.
 */
//---------------------------------------------------------------------------//

#include "DTK_STKMeshManager.hpp"
#include "DTK_STKMeshEntityPredicates.hpp"
#include "DTK_STKMeshNodalShapeFunction.hpp"
#include "DTK_STKMeshEntitySet.hpp"
#include "DTK_STKMeshEntityLocalMap.hpp"
#include "DTK_STKMeshEntityIntegrationRule.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
//! Default constructor.
STKMeshManager::STKMeshManager(
    const Teuchos::RCP<stk::mesh::BulkData>& bulk_data,
    const BasisType basis_type )
    : d_bulk_data( bulk_data )
{
    createFunctionSpace( basis_type, FunctionSpace::selectAll );
    DTK_ENSURE( Teuchos::nonnull(d_function_space) );
}

//---------------------------------------------------------------------------//
//! Part name constructor.
STKMeshManager::STKMeshManager(
    const Teuchos::RCP<stk::mesh::BulkData>& bulk_data,
    const Teuchos::Array<std::string>& part_names,
    const BasisType basis_type )
    : d_bulk_data( bulk_data )
{
    STKPartNamePredicate pred( part_names, d_bulk_data );
    createFunctionSpace( basis_type, pred.getFunction() );
    DTK_ENSURE( Teuchos::nonnull(d_function_space) );
}

//---------------------------------------------------------------------------//
//! Part vector constructor.
STKMeshManager::STKMeshManager(
    const Teuchos::RCP<stk::mesh::BulkData>& bulk_data,
    const stk::mesh::PartVector& parts,
    const BasisType basis_type )
    : d_bulk_data( bulk_data )
{
    STKPartVectorPredicate pred( parts );
    createFunctionSpace( basis_type, pred.getFunction() );
    DTK_ENSURE( Teuchos::nonnull(d_function_space) );
}

//---------------------------------------------------------------------------//
//! Selector constructor.
STKMeshManager::STKMeshManager(
    const Teuchos::RCP<stk::mesh::BulkData>& bulk_data,
    const stk::mesh::Selector& selector,
    const BasisType basis_type )
    : d_bulk_data( bulk_data )
{
    STKSelectorPredicate pred( selector );
    createFunctionSpace( basis_type, pred.getFunction() );
    DTK_ENSURE( Teuchos::nonnull(d_function_space) );
}

//---------------------------------------------------------------------------//
// Get the function space over which the mesh and its fields are defined.
Teuchos::RCP<FunctionSpace> STKMeshManager::functionSpace() const
{
    return d_function_space;
}

//---------------------------------------------------------------------------//
// Create the function space.
void STKMeshManager::createFunctionSpace(
    const BasisType basis_type,
    const PredicateFunction& select_function )
{
    Teuchos::RCP<EntitySet> entity_set =
        Teuchos::rcp( new STKMeshEntitySet(d_bulk_data) );

    Teuchos::RCP<EntityLocalMap> local_map =
        Teuchos::rcp( new STKMeshEntityLocalMap(d_bulk_data) );

    Teuchos::RCP<EntityShapeFunction> shape_function;
    switch( basis_type )
    {
        case BASIS_TYPE_GRADIENT:
            shape_function =
                Teuchos::rcp( new STKMeshNodalShapeFunction(d_bulk_data) );
            break;

        default:
            bool bad_basis_type = true;
            DTK_INSIST( !bad_basis_type );
            break;
    }
    DTK_CHECK( Teuchos::nonnull(shape_function) );

    Teuchos::RCP<EntityIntegrationRule> integration_rule =
        Teuchos::rcp( new STKMeshEntityIntegrationRule(d_bulk_data) );

    d_function_space = Teuchos::rcp(
        new FunctionSpace(entity_set,local_map,shape_function,
                          integration_rule,select_function) );

    DTK_ENSURE( Teuchos::nonnull(d_function_space) );
}

//---------------------------------------------------------------------------//
Teuchos::RCP<EntitySet> STKMeshManager::entitySet() const
{
    return d_function_space->entitySet();
}

//---------------------------------------------------------------------------//
// Get the local map for entities supporting the function.
Teuchos::RCP<EntityLocalMap> STKMeshManager::localMap() const
{
    return d_function_space->localMap();
}

//---------------------------------------------------------------------------//
// Get the shape function for entities supporting the function.
Teuchos::RCP<EntityShapeFunction> STKMeshManager::shapeFunction() const
{
    return d_function_space->shapeFunction();
}

//---------------------------------------------------------------------------//
// Get the integration rule for entities supporting the function.
Teuchos::RCP<EntityIntegrationRule> STKMeshManager::integrationRule() const
{
    return d_function_space->integrationRule();
}

//---------------------------------------------------------------------------//
// Get the selector function.
PredicateFunction STKMeshManager::selectFunction() const
{
    return d_function_space->selectFunction();
}

//---------------------------------------------------------------------------//
// Get the field for the given string key.
Teuchos::RCP<Field>
STKMeshManager::field( const std::string& field_name ) const
{
    DTK_REQUIRE( d_field_indexer.count(field_name) );
    int field_id = d_field_indexer.find(field_name)->second;
    return d_fields[ field_id ];
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_STKMeshManager.cpp
//---------------------------------------------------------------------------//
