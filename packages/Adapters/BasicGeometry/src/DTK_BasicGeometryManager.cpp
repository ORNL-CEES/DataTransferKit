//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \brief DTK_BasicGeometryManager.cpp
 * \author Stuart R. Slattery
 * \brief High-level manager for basic geometries.
 */
//---------------------------------------------------------------------------//

#include "DTK_BasicEntityPredicates.hpp"
#include "DTK_BasicGeometryManager.hpp"
#include "DTK_BasicEntitySet.hpp"
#include "DTK_BasicGeometryLocalMap.hpp"
#include "DTK_EntityCenteredShapeFunction.hpp"
#include "DTK_PredicateComposition.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
//! Default constructor.
BasicGeometryManager::BasicGeometryManager(
    const Teuchos::RCP<const Teuchos::Comm<int> > comm,
    const int physical_dimension,
    const EntityType entity_type )
{
    Teuchos::RCP<EntitySelector> entity_selector = 
	Teuchos::rcp( new EntitySelector(entity_type) );
    Teuchos::Array<Entity> entities(0);
    createFunctionSpace( comm, physical_dimension, entities, entity_selector );
    DTK_ENSURE( Teuchos::nonnull(d_function_space) );
}

//---------------------------------------------------------------------------//
//! Entity constructor.
BasicGeometryManager::BasicGeometryManager(
    const Teuchos::RCP<const Teuchos::Comm<int> > comm,
    const int physical_dimension,
    const EntityType entity_type,
    const Teuchos::ArrayView<Entity>& entities )
{
    Teuchos::RCP<EntitySelector> entity_selector = 
	Teuchos::rcp( new EntitySelector(entity_type) );
    createFunctionSpace( comm, physical_dimension, entities, entity_selector );
    DTK_ENSURE( Teuchos::nonnull(d_function_space) );
}

//---------------------------------------------------------------------------//
//! Predicate constructor.
BasicGeometryManager::BasicGeometryManager(
    const Teuchos::RCP<const Teuchos::Comm<int> > comm,
    const int physical_dimension,
    const EntityType entity_type,
    const Teuchos::ArrayView<Entity>& entities,
    const Teuchos::ArrayView<int>& block_ids,
    const Teuchos::ArrayView<int>& boundary_ids )
{
    Teuchos::Array<int> block_array(block_ids);
    Teuchos::Array<int> boundary_array(boundary_ids);
    BlockPredicate block_pred( block_array );
    BoundaryPredicate boundary_pred( boundary_array );
    std::function<bool(Entity)> pred = 
	PredicateComposition::Or( block_pred.getFunction(),
				  boundary_pred.getFunction() );
    Teuchos::RCP<EntitySelector> entity_selector = 
	Teuchos::rcp( new EntitySelector(entity_type,pred) );
    createFunctionSpace( comm, physical_dimension, entities, entity_selector );
    DTK_ENSURE( Teuchos::nonnull(d_function_space) );
}

//---------------------------------------------------------------------------//
// Get the function space over which the mesh and its fields are defined. 
Teuchos::RCP<FunctionSpace> BasicGeometryManager::functionSpace() const
{
    return d_function_space;
}

//---------------------------------------------------------------------------//
// Create the function space.
void BasicGeometryManager::createFunctionSpace( 
	const Teuchos::RCP<const Teuchos::Comm<int> > comm,
	const int physical_dimension, 
	const Teuchos::ArrayView<Entity>& entities,
	const Teuchos::RCP<EntitySelector>& entity_selector )
{
    Teuchos::RCP<BasicEntitySet> entity_set = 
	Teuchos::rcp( new BasicEntitySet(comm,physical_dimension) );
    for ( Entity e : entities )
    {
	entity_set->addEntity( e );
    }
    
    Teuchos::RCP<EntityLocalMap> local_map =
	Teuchos::rcp( new BasicGeometryLocalMap() );

    Teuchos::RCP<EntityShapeFunction> shape_function =
	Teuchos::rcp( new EntityCenteredShapeFunction() );

    d_function_space = Teuchos::rcp( 
	new FunctionSpace(entity_set,entity_selector,local_map,shape_function) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_BasicGeometryManager.cpp
//---------------------------------------------------------------------------//
