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
 * \brief DTK_MoabManager.cpp
 * \author Stuart R. Slattery
 * \brief High-level manager for Moab.
 */
//---------------------------------------------------------------------------//

#include "DTK_MoabManager.hpp"
#include "DTK_MoabEntityPredicates.hpp"
#include "DTK_MoabNodalShapeFunction.hpp"
#include "DTK_MoabEntitySet.hpp"
#include "DTK_MoabEntityLocalMap.hpp"
#include "DTK_MoabMeshSetIndexer.hpp"
#include "DTK_MoabEntityIntegrationRule.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
//! Default constructor.
MoabManager::MoabManager( const Teuchos::RCP<moab::ParallelComm>& moab_mesh,
			  bool create_global_ids )
    : d_moab_mesh( moab_mesh )
{
    // Build a mesh set indexer. This must be constructed after the global
    // indices are created.
    d_set_indexer =
	Teuchos::rcp( new MoabMeshSetIndexer(d_moab_mesh,create_global_ids) );
    
    // Build DTK data structures.
    Teuchos::RCP<EntitySet> entity_set = 
	Teuchos::rcp( new MoabEntitySet(d_moab_mesh,d_set_indexer) );
    
    Teuchos::RCP<EntityLocalMap> local_map =
	Teuchos::rcp( new MoabEntityLocalMap(d_moab_mesh) );

    Teuchos::RCP<EntityShapeFunction> shape_function =
	Teuchos::rcp( new MoabNodalShapeFunction(d_moab_mesh) );

    Teuchos::RCP<EntityIntegrationRule> integration_rule =
	Teuchos::rcp( new MoabEntityIntegrationRule(d_moab_mesh) );

    d_function_space = Teuchos::rcp( 
	new FunctionSpace(entity_set,local_map,shape_function,integration_rule) );
}

//---------------------------------------------------------------------------//
//! Mesh set constructor.
MoabManager::MoabManager( const Teuchos::RCP<moab::ParallelComm>& moab_mesh,
			  const moab::EntityHandle& mesh_set,
			  bool create_global_ids )
    : d_moab_mesh( moab_mesh )
{
    // Build a mesh set indexer. This must be constructed after the global
    // indices are created.
    d_set_indexer =
	Teuchos::rcp( new MoabMeshSetIndexer(d_moab_mesh,create_global_ids) );

    // Build DTK data structures.
    Teuchos::RCP<EntitySet> entity_set = 
	Teuchos::rcp( new MoabEntitySet(d_moab_mesh,d_set_indexer) );

    MoabMeshSetPredicate pred( mesh_set, d_set_indexer );

   Teuchos::RCP<EntityLocalMap> local_map =
	Teuchos::rcp( new MoabEntityLocalMap(d_moab_mesh) );

    Teuchos::RCP<EntityShapeFunction> shape_function =
	Teuchos::rcp( new MoabNodalShapeFunction(d_moab_mesh) );

    Teuchos::RCP<EntityIntegrationRule> integration_rule =
	Teuchos::rcp( new MoabEntityIntegrationRule(d_moab_mesh) );

    d_function_space = Teuchos::rcp( new FunctionSpace(entity_set,
						       local_map,
						       shape_function,
						       integration_rule,
						       pred.getFunction()) );
}

//---------------------------------------------------------------------------//
// Register a tag and associated entity set with the manager that will be
// available for solution transfer. 
void MoabManager::registerTag( const moab::EntityHandle& mesh_set,
			       const moab::Tag& tag )
{
    std::string tag_name;
    DTK_CHECK_ERROR_CODE(
	d_moab_mesh->get_moab()->tag_get_name( tag, tag_name )
	);
    d_field_indexer.emplace( tag_name, d_fields.size() );
    d_fields.push_back(
	Teuchos::rcp( new MoabTagField<double>(d_moab_mesh,
					       d_set_indexer,
					       mesh_set,
					       tag) )
	);
}

//---------------------------------------------------------------------------//
// Get the function space over which the mesh and its fields are defined. 
Teuchos::RCP<FunctionSpace> MoabManager::functionSpace() const
{
    return d_function_space;
}

//---------------------------------------------------------------------------//
Teuchos::RCP<FieldMultiVector<double> >
MoabManager::createFieldMultiVector( const moab::EntityHandle& mesh_set,
				     const moab::Tag& tag )
{
    DTK_REQUIRE( Teuchos::nonnull(d_moab_mesh) );
    DTK_REQUIRE( Teuchos::nonnull(d_function_space) );
    
    Teuchos::RCP<Field<double> > field = Teuchos::rcp(
	new MoabTagField<double>(d_moab_mesh, d_set_indexer, mesh_set, tag) );
    return Teuchos::rcp(
	new FieldMultiVector<double>(field,d_function_space->entitySet()) );
}

//---------------------------------------------------------------------------//
Teuchos::RCP<EntitySet> MoabManager::entitySet() const
{
    return d_function_space->entitySet();
}

//---------------------------------------------------------------------------//
// Get the local map for entities supporting the function.
Teuchos::RCP<EntityLocalMap> MoabManager::localMap() const
{
    return d_function_space->localMap();
}

//---------------------------------------------------------------------------//
// Get the shape function for entities supporting the function.
Teuchos::RCP<EntityShapeFunction> MoabManager::shapeFunction() const
{
    return d_function_space->shapeFunction();
}

//---------------------------------------------------------------------------//
// Get the integration rule for entities supporting the function.
Teuchos::RCP<EntityIntegrationRule> MoabManager::integrationRule() const
{
    return d_function_space->integrationRule();
}

//---------------------------------------------------------------------------//
// Get the selector function.
PredicateFunction MoabManager::selectFunction() const
{
    return d_function_space->selectFunction();
}

//---------------------------------------------------------------------------//
// Get the field for the given string key.
Teuchos::RCP<Field<double> >
MoabManager::field( const std::string& field_name ) const
{
    DTK_REQUIRE( d_field_indexer.count(field_name) );
    int field_id = d_field_indexer.find(field_name)->second;
    return d_fields[ field_id ];
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MoabManager.cpp
//---------------------------------------------------------------------------//
