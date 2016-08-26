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
 * \brief DTK_IntrepidShapeFunction.cpp
 * \author Stuart R. Slattery
 * \brief  Intrepid shape function implementation.
 */
//---------------------------------------------------------------------------//

#include "DTK_IntrepidShapeFunction.hpp"
#include "DTK_IntrepidBasisFactory.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Given an topology and a reference point, evaluate the shape function of the
// topology at that point.
void IntrepidShapeFunction::evaluateValue( 
    const shards::CellTopology& topology,
    const Teuchos::ArrayView<const double>& reference_point,
    Teuchos::Array<double>& values ) const
{
    // Get the basis for the topology.
    Teuchos::RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
	basis = getIntrepidBasis( topology );

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
// Given an topology and a reference point, evaluate the gradient of the shape
// function of the topology at that point.
void IntrepidShapeFunction::evaluateGradient( 
    const shards::CellTopology& topology,
    const Teuchos::ArrayView<const double>& reference_point,
    Teuchos::Array<Teuchos::Array<double> >& gradients ) const
{
    // Get the basis for the topology.
    Teuchos::RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
	basis = getIntrepidBasis( topology );

    // Wrap the reference point.
    int space_dim = reference_point.size();
    Teuchos::Array<int> point_dims(2);
    point_dims[0] = 1;
    point_dims[1] = space_dim;
    Intrepid::FieldContainer<double> point_container( 
	point_dims, const_cast<double*>(reference_point.getRawPtr()) );

    // Evaluate the basis function.
    int cardinality = basis->getCardinality();
    Intrepid::FieldContainer<double> grad_container( cardinality, 1, space_dim );
    basis->getValues( 
	grad_container, point_container, Intrepid::OPERATOR_GRAD );

    // Extract the evaluations.
    gradients.resize( cardinality );
    for ( int n = 0; n < cardinality; ++n )
    {
	gradients[n].resize( space_dim );
	for ( int d = 0; d < space_dim; ++d )
	{
	    gradients[n][d] = grad_container(n,0,d);
	}
    }
}

//---------------------------------------------------------------------------//
// Given a topology, get the intrepid basis function.
Teuchos::RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
IntrepidShapeFunction::getIntrepidBasis(
    const shards::CellTopology& topology ) const
{
    // Either make a new basis for this topology or return an existing one.
    unsigned basis_key = topology.getKey();
    Teuchos::RCP<
        Intrepid::Basis<double,Intrepid::FieldContainer<double> > > basis;
    if ( d_basis.count(basis_key) )
    {
	basis = d_basis.find( basis_key )->second;
    }
    else
    {
	basis = IntrepidBasisFactory::create( topology );
	d_basis.emplace( basis_key, basis );
    }
    return basis;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_IntrepidShapeFunction.cpp
//---------------------------------------------------------------------------//
