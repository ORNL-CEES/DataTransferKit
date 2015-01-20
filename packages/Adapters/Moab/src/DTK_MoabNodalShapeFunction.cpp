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
 * \brief DTK_MoabNodalShapeFunction.cpp
 * \author Stuart R. Slattery
 * \brief Nodal shape function implementation for STK mesh.
 */
//---------------------------------------------------------------------------//

#include "DTK_MoabNodalShapeFunction.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
MoabNodalShapeFunction::MoabNodalShapeFunction(
    const Teuchos::RCP<moab::ParallelComm>& moab_mesh )
    : d_moab_mesh( moab_mesh )
{
    d_moab_evaluator = Teuchos::rcp( 
	new moab::ElemEvaluator(d_moab_mesh->get_moab()) );
}

//---------------------------------------------------------------------------//
// Destructor.
MoabNodalShapeFunction::~MoabNodalShapeFunction()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Given an entity, get the ids of the degrees of freedom in the vector space
// supporting its shape function.
void MoabNodalShapeFunction::entityDOFIds( 
    const Entity& entity, Teuchos::Array<std::size_t>& dof_ids ) const
{
    // Node case.
    if ( ENTITY_TYPE_NODE == entity.entityType() )
    {
	dof_ids.assign( 1, entity.id() );
    }

    // Element case.
    else
    {
	const moab::EntityHandle* entity_nodes;
	int num_nodes = 0;
	std::vector<moab::EntityHandle> storage;
	DTK_CHECK_ERROR_CODE(
	    d_moab_mesh->get_moab()->get_connectivity( entity.id(),
						       entity_nodes,
						       num_nodes,
						       false,
						       &storage )
	    );

	dof_ids.resize( num_nodes );
	for ( int i = 0; i < num_nodes; ++i )
	{
	    dof_ids[i] = Teuchos::as<std::size_t>( entity_nodes[i] );
	}
    }
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the shape function of the
// entity at that point.
void MoabNodalShapeFunction::evaluateValue( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    Teuchos::Array<double>& values ) const
{
    // Cache the entity with the evaluator.
    cacheEntity( entity );

    // Get the number of nodes supporting the entity.
    const moab::EntityHandle* entity_nodes;
    int num_nodes = 0;
    DTK_CHECK_ERROR_CODE(
	d_moab_mesh->get_moab()->get_connectivity( entity.id(),
						   entity_nodes,
						   num_nodes )
	);

    values.resize( num_nodes );

    // Extract the value of the basis used by the evaluator by passing the
    // identity matrix through the eval function.
    int topo_dim =
	d_moab_mesh->get_moab()->dimension_from_handle( entity.id() );
    Teuchos::Array<double> field( num_nodes*num_nodes, 0.0 );
    for ( int n = 0; n < num_nodes; ++n )
    {
	field[n*num_nodes + n] = 1.0;
    }
 
    moab::EntityType moab_type = 
	d_moab_mesh->get_moab()->type_from_handle( entity.id() );
    DTK_CHECK_ERROR_CODE(
	(*d_moab_evaluator->get_eval_set(moab_type).evalFcn)
	( reference_point.getRawPtr(),
	  field.getRawPtr(),
	  topo_dim,
	  num_nodes,
	  d_moab_evaluator->get_work_space(),
	  values.getRawPtr() )
	);
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the gradient of the shape
// function of the entity at that point.
void MoabNodalShapeFunction::evaluateGradient( 
	const Entity& entity,
	const Teuchos::ArrayView<const double>& reference_point,
	Teuchos::Array<Teuchos::Array<double> >& gradients ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Cache an entity in the evaluator.
void MoabNodalShapeFunction::cacheEntity( const Entity& entity ) const
{
    DTK_CHECK_ERROR_CODE(
	d_moab_evaluator->set_eval_set( entity.id() )
	);
    DTK_CHECK_ERROR_CODE(
	d_moab_evaluator->set_ent_handle( entity.id() )
	);
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MoabNodalShapeFunction.cpp
//---------------------------------------------------------------------------//
