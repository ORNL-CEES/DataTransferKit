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
 * \file   DTK_LocalMLSProblem_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Local moving least square problem.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_LOCALMLSPROBLEM_IMPL_HPP
#define DTK_LOCALMLSPROBLEM_IMPL_HPP

#include "DTK_RadialBasisPolicy.hpp"
#include "DTK_EuclideanDistance.hpp"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Basis, class GO, int DIM>
LocalMLSProblem<Basis,GO,DIM>::LocalMLSProblem( 
    const Teuchos::ArrayView<const double>& target_center,
    const Teuchos::ArrayView<const unsigned>& source_lids,
    const Teuchos::ArrayView<const double>& source_centers,
    const Basis& basis,
    const double alpha )
    : d_shape_function( source_lids.size() )
{
    DTK_REQUIRE( 0 == source_centers.size() % DIM );
    DTK_REQUIRE( 0 == target_center.size() % DIM );
    DTK_REQUIRE( alpha >= 0.0 );

    // Number of source centers supporting this target center.
    int num_sources = source_lids.size();

    // Build the matrix of basis evaluations and the P matrix.
    int poly_size = 0;
    if ( 1 == DIM ) poly_size = 3;
    else if ( 2 == DIM ) poly_size = 6;
    else if ( 3 == DIM ) poly_size = 10;
    Teuchos::SerialDenseMatrix<int,double> P( num_sources, poly_size );
    Teuchos::SerialDenseMatrix<int,double> phi( num_sources, num_sources );
    Teuchos::ArrayView<const double> source_center_view;
    double dist = 0.0;
    for ( int i = 0; i < num_sources; ++i )
    {
	source_center_view = source_centers(DIM*source_lids[i],DIM);
	dist = EuclideanDistance<DIM>::distance(
	    target_center.getRawPtr(), source_center_view.getRawPtr() );

	// Basis values.
	phi(i,i) = BP::evaluateValue( basis, dist );
	if ( alpha > 0.0 )
	{
	    phi(i,i) += alpha * BP::evaluateGradient( basis, dist );
	}

	// Polynomial matrix.
	P(i,0) = 1;
	P(i,1) = source_center_view[0];
	if ( 1 == DIM )
	{
	    P(i,2) = source_center_view[0]*source_center_view[0];
	}
	else if ( 2 == DIM )
	{
	    P(i,2) = source_center_view[1];
	    P(i,3) = source_center_view[0]*source_center_view[0];
	    P(i,4) = source_center_view[0]*source_center_view[1];
	    P(i,5) = source_center_view[1]*source_center_view[1];
	}
	else if ( 3 == DIM )
	{
	    P(i,2) = source_center_view[1];
	    P(i,3) = source_center_view[2];
	    P(i,4) = source_center_view[0]*source_center_view[0];
	    P(i,5) = source_center_view[0]*source_center_view[1];
	    P(i,6) = source_center_view[1]*source_center_view[1];
	    P(i,7) = source_center_view[1]*source_center_view[2];
	    P(i,8) = source_center_view[2]*source_center_view[2];
	    P(i,9) = source_center_view[2]*source_center_view[0];
	}
    }

    // Build the target polynomial.
    Teuchos::SerialDenseVector<int,double> target_poly( poly_size );
    target_poly(0) = 1;
    target_poly(1) = target_center[0];
    if ( 1 == DIM )
    {
	target_poly(2) = target_center[0]*target_center[0];
    }
    else if ( 2 == DIM )
    {
	target_poly(2) = target_center[1];
	target_poly(3) = target_center[0]*target_center[0];
	target_poly(4) = target_center[0]*target_center[1];
	target_poly(5) = target_center[1]*target_center[1];
    }
    else if ( 3 == DIM )
    {
	target_poly(2) = target_center[1];
	target_poly(3) = target_center[2];
	target_poly(4) = target_center[0]*target_center[0];
	target_poly(5) = target_center[0]*target_center[1];
	target_poly(6) = target_center[1]*target_center[1];
	target_poly(7) = target_center[1]*target_center[2];
	target_poly(8) = target_center[2]*target_center[2];
	target_poly(9) = target_center[2]*target_center[0];
    }

    // Construct and invert the A matrix.
    Teuchos::SerialDenseMatrix<int,double> A( poly_size, poly_size );
    {
	Teuchos::SerialDenseMatrix<int,double> work( num_sources, poly_size );
	work.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, phi, P, 0.0 );
	A.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, P, work, 0.0 );
	Teuchos::SerialDenseSolver<int,double> A_solver;
	A_solver.setMatrix( Teuchos::rcpFromRef(A) );
	A_solver.invert();
    }

    // Construct the basis.
    Teuchos::SerialDenseMatrix<int,double> b( poly_size, num_sources );
    b.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, P, phi, 0.0 );
    Teuchos::SerialDenseMatrix<int,double> work( poly_size, num_sources );
    work.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, b, 0.0 );
    Teuchos::SerialDenseMatrix<int,double> shape_matrix(
	Teuchos::View, d_shape_function.getRawPtr(), 
	1, 1, d_shape_function.size() );
    shape_matrix.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 
			   1.0, target_poly, work, 0.0 );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the local problem by applying the shape function to the
 * degrees of freedom at the given source centers.
 */
template<class Basis, class GO, int DIM>
double LocalMLSProblem<Basis,GO,DIM>::solve( 
    const Teuchos::ArrayView<const double>& source_data,
    const Teuchos::ArrayView<const unsigned>& source_lids ) const
{
    DTK_REQUIRE( source_lids.size() == d_shape_function.size() );
    double target_data = 0.0;
    Teuchos::Array<double>::const_iterator shape_it;
    Teuchos::ArrayView<const unsigned>::const_iterator lid_it;
    for ( shape_it = d_shape_function.begin(), lid_it = source_lids.begin();
	  shape_it != d_shape_function.end();
	  ++shape_it, ++lid_it )
    {
	target_data += *shape_it * source_data[*lid_it];
    }

    return target_data;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_LOCALMLSPROBLEM_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_LocalMLSProblem_impl.hpp
//---------------------------------------------------------------------------//

