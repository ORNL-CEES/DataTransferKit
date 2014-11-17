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
 * \brief DTK_STKMeshNodalShapeFunction.cpp
 * \author Stuart R. Slattery
 * \brief Nodal shape function implementation for STK mesh.
 */
//---------------------------------------------------------------------------//

#include "DTK_STKMeshNodalShapeFunction.hpp"
#include "DTK_STKMeshHelpers.hpp"
#include "DTK_IntrepidBasisFactory.hpp"
#include "DTK_DBC.hpp"

#include <Shards_CellTopology.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
STKMeshNodalShapeFunction::STKMeshNodalShapeFunction(
    const Teuchos::RCP<stk::mesh::BulkData>& bulk_data )
    : d_bulk_data( bulk_data )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
STKMeshNodalShapeFunction::~STKMeshNodalShapeFunction()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Given an entity, get the ids of the degrees of freedom in the vector space
// supporting its shape function.
void STKMeshNodalShapeFunction::entityDOFIds( 
    const Entity& entity, Teuchos::Array<std::size_t>& dof_ids ) const
{
    // Extract the stk entity.
    const stk::mesh::Entity& stk_entity = 
	STKMeshHelpers::extractEntity( entity );

    // Get the ids of the nodes supporting the entity.
    const stk::mesh::Entity* begin = d_bulk_data->begin_nodes( stk_entity );
    const stk::mesh::Entity* end = d_bulk_data->end_nodes( stk_entity );

    // Extract the node ids as the dof ids.
    int num_nodes = std::distance( begin, end );
    dof_ids.resize( num_nodes );
    for ( int n = 0; n < num_nodes; ++n )
    {
	dof_ids[n] = d_bulk_data->identifier( begin[n] );
    }
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the shape function of the
// entity at that point.
void STKMeshNodalShapeFunction::evaluateValue( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    Teuchos::Array<double>& values ) const
{
    // Get the basis for the entity.
    Teuchos::RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
	basis = getIntrepidBasis( entity );

    // Wrap the reference point.
    Teuchos::Array<int> point_dims(2);
    point_dims[0] = 1;
    point_dims[1] = reference_point.size();
    Intrepid::FieldContainer<double> point_container(
	point_dims, const_cast<double*>(reference_point.getRawPtr()) );

    // Wrap the evaluations.
    values.resize( basis->getCardinality() );
    Teuchos::Array<int> value_dims(2);
    value_dims[0] = basis->getCardinality();
    value_dims[1] = 1;
    Intrepid::FieldContainer<double> value_container(
	value_dims, values.getRawPtr() );

    // Evaluate the basis function.
    basis->getValues( 
	value_container, point_container, Intrepid::OPERATOR_VALUE );
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the gradient of the shape
// function of the entity at that point.
void STKMeshNodalShapeFunction::evaluateGradient( 
	const Entity& entity,
	const Teuchos::ArrayView<const double>& reference_point,
	Teuchos::Array<Teuchos::Array<double> >& gradients ) const
{
    // Get the basis for the entity.
    Teuchos::RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
	basis = getIntrepidBasis( entity );

    // Wrap the reference point.
    int space_dim = reference_point.size();
    Teuchos::Array<int> point_dims(2);
    point_dims[0] = 1;
    point_dims[1] = space_dim;
    Intrepid::FieldContainer<double> point_container( 
	point_dims, const_cast<double*>(reference_point.getRawPtr()) );

    // Evaluate the basis function.
    int cardinality = basis->getCardinality();
    Intrepid::FieldContainer<double> grad_container( cardinality, space_dim );
    basis->getValues( 
	grad_container, point_container, Intrepid::OPERATOR_GRAD );

    // Extract the evaluations.
    gradients.resize( cardinality );
    for ( int n = 0; n < cardinality; ++n )
    {
	gradients[n].resize( space_dim );
	for ( int d = 0; d < space_dim; ++d )
	{
	    gradients[n][d] = grad_container(n,d);
	}
    }
}

//---------------------------------------------------------------------------//
// Given an entity, get the intrepid basis function.
Teuchos::RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
STKMeshNodalShapeFunction::getIntrepidBasis( const Entity& entity ) const
{
    // Extract the stk entity.
    const stk::mesh::Entity& stk_entity = 
	STKMeshHelpers::extractEntity( entity );

    // Get the topology of the entity.
    shards::CellTopology entity_topo = 
	stk::mesh::get_cell_topology(
	    d_bulk_data->bucket(stk_entity).topology() );

    // Get the basis.
    return IntrepidBasisFactory::create( entity_topo );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_STKMeshNodalShapeFunction.cpp
//---------------------------------------------------------------------------//
