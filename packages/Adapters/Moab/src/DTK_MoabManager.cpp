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

    d_function_space = Teuchos::rcp( 
	new FunctionSpace(entity_set,local_map,shape_function) );
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

    d_function_space = Teuchos::rcp(
	new FunctionSpace(entity_set,
			  local_map,
			  shape_function,
			  pred.getFunction()) );
}

//---------------------------------------------------------------------------//
// Get the function space over which the mesh and its fields are defined. 
Teuchos::RCP<FunctionSpace> MoabManager::functionSpace() const
{
    return d_function_space;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MoabManager.cpp
//---------------------------------------------------------------------------//
