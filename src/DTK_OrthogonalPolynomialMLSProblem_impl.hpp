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
 * \file   DTK_OrthogonalPolynomialMLSProblem_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Local moving least square problem.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ORTHOGONALPOLYNOMIALMLSPROBLEM_IMPL_HPP
#define DTK_ORTHOGONALPOLYNOMIALMLSPROBLEM_IMPL_HPP

#include <limits>

#include "DTK_DBC.hpp"
#include "DTK_RadialBasisPolicy.hpp"
#include "DTK_EuclideanDistance.hpp"

#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_LAPACK.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Basis, class GO, int DIM>
OrthogonalPolynomialMLSProblem<Basis,GO,DIM>::OrthogonalPolynomialMLSProblem( 
    const Teuchos::ArrayView<const double>& target_center,
    const Teuchos::ArrayView<const unsigned>& source_lids,
    const Teuchos::ArrayView<const double>& source_centers,
    const Basis& basis )
    : d_shape_function( source_lids.size() )
{
    DTK_REQUIRE( 0 == source_centers.size() % DIM );
    DTK_REQUIRE( 0 == target_center.size() % DIM );

    // Number of source centers supporting this target center.
    int num_sources = source_lids.size();
    DTK_CHECK( 0 < num_sources );

    // Build the matrix of basis evaluations.
    Teuchos::Array<double> Phi( num_sources ); 
    double dist = 0.0;
    Teuchos::ArrayView<const double> source_center_view;
    for ( int i = 0; i < num_sources; ++i )
    {
	source_center_view = source_centers(DIM*source_lids[i],DIM);
	dist = EuclideanDistance<DIM>::distance(
	    target_center.getRawPtr(), source_center_view.getRawPtr() );
	Phi[i] = BP::evaluateValue( basis, dist );
    }

    // Build the polynomial space of the supporting source points.
    Teuchos::Array<Teuchos::Array<double> > source_space( num_sources );
    for ( int i = 0; i < num_sources; ++i )
    {
	source_space[i] = 
	    polynomialSpace( source_centers(DIM*source_lids[i],DIM) );
    }
    int space_size = source_space[0].size();
    DTK_CHECK( space_size >= num_sources );

    // Build the matrix of monomials.
    Teuchos::Array<int> space_components( num_sources, -1 );
    Teuchos::SerialDenseMatrix<int,double> X;
    Teuchos::LAPACK<int,double> lapack;
    int info = 0;
    {
	Teuchos::SerialDenseMatrix<int,double> v( num_sources, 1 );
	Teuchos::SerialDenseVector<int,double> s;
	Teuchos::Array<double> work( 4 * num_sources );
	int rank = 0;
	int space_id = 0;
	bool is_full_rank = false;
	for ( int i = 0; i < num_sources; ++i )
	{
	    // Add a row to X.
	    X.reshape( i+1, num_sources );

	    // Get the optimal work size for GELSS.
	    s.resize( i+1 );
	    Teuchos::SerialDenseMatrix<int,double> Xcp( X );
	    lapack.GELSS( Xcp.numRows(), Xcp.numCols(), v.numCols(), 
			  Xcp.values(), Xcp.numRows(),
			  v.values(), v.numRows(), s.values(),
			  10.0*std::numeric_limits<double>::epsilon(), 
			  &rank, work.getRawPtr(), -1, &info );
	    work.resize( work[0] );

	    // Cycle through the polynomial space until we end up with a new X
	    // that is not rank deficient or we run out of polynomial
	    // components. If we run out this will trigger a DBC check that is
	    // always on regardless of DBC build status.
	    is_full_rank = false;
	    while ( !is_full_rank && space_id < space_size )
	    {
		// Assign the next polynomial space to the last row.
		DTK_CHECK( space_id < space_size );
		for ( int j = 0; j < num_sources; ++j )
		{
		    X(i,j) = source_space[j][space_id];
		}

		// Check if X is full rank using SVD.
		Xcp = Teuchos::SerialDenseMatrix<int,double>( X );
		lapack.GELSS( Xcp.numRows(), Xcp.numCols(), v.numCols(), 
			      Xcp.values(), Xcp.numRows(),
			      v.values(), v.numRows(), s.values(),
			      -1.0, &rank, work.getRawPtr(), work.size(), &info );
		DTK_CHECK( 0 == info );
		is_full_rank = ( rank == i+1 );
		std::cout << s << std::endl;
		// Update the polynomial space component id.
		++space_id;
	    }
	    DTK_CHECK( is_full_rank );

	    // Store the polynomial space component that made this row.
	    space_components[i] = space_id-1;
	}
	DTK_REQUIRE( is_full_rank );
	DTK_CHECK( space_id < space_size );
	DTK_CHECK( X.numRows() == num_sources );
	DTK_CHECK( X.numCols() == num_sources );
    }

    // Build the M matrix.
    Teuchos::SerialDenseMatrix<int,double> M( num_sources, num_sources );
    M.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, X, X, 0.0 );

    // Compute the LDLT factorization of M.
    factorLDLT( M );

    // Extract the diagonal and scale.
    Teuchos::SerialDenseMatrix<int,double> D( num_sources, num_sources );
    for ( int i = 0; i < num_sources; ++i )
    {
	D(i,i) = 1.0 / std::sqrt( M(i,i)*Phi[i] );
    }

    // Extract the lower triangular factor transpose L^T.
    Teuchos::SerialDenseMatrix<int,double> LT( num_sources, num_sources );
    for ( int i = 0; i < num_sources; ++i )
    {
	LT(i,i) = 1.0;
	for ( int j = 0; j < num_sources; ++j )
	{
	    if ( i > j )
	    {
		LT(j,i) = M(i,j);
	    }
	}
    }

    // Compute the polynomial coefficients via a triangular solve.
    lapack.TRTRS( 'U', 'N', 'U', LT.numRows(), D.numCols(), LT.values(),
		  LT.numRows(), D.values(), D.numRows(), &info );
    DTK_CHECK( 0 == info );

    // Build the target polynomial.
    Teuchos::Array<double> target_space = polynomialSpace( target_center );
    Teuchos::SerialDenseVector<int,double> target_poly( num_sources );
    for ( int i = 0; i < num_sources; ++i )
    {
	target_poly(i) = target_space[ space_components[i] ];
    }
    Teuchos::SerialDenseMatrix<int,double> Pt( num_sources, 1 );
    Pt.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, D, target_poly, 0.0 );

    // Build the source polynomial.
    Teuchos::SerialDenseMatrix<int,double> source_poly( num_sources, num_sources );
    source_poly.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, D, X, 0.0 );
    std::cout << D << std::endl;

    // Construct the basis.
    for ( int j = 0; j < num_sources; ++j )
    {
	for ( int i = 0; i < num_sources; ++i )
	{
	    d_shape_function[j] += source_poly(i,j)*target_poly(i);
	}
	d_shape_function[j] *= Phi[j];
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given a set of coordinates, compute a polynomial space.
 */
template<class Basis, class GO, int DIM>
Teuchos::Array<double> OrthogonalPolynomialMLSProblem<Basis,GO,DIM>::polynomialSpace(
    const Teuchos::ArrayView<const double>& center ) const
{
    Teuchos::Array<double> poly_space;

    // 1 dimension.
    if ( 1 == DIM )
    {
	int poly_size = 20;
	poly_space.resize( poly_size );
	for ( int i = 0; i < poly_size; ++i )
	{
	    poly_space[i] = std::pow( center[0], i );
	}
    }

    // 2 dimensions.
    else if ( 2 == DIM )
    {
	int poly_size = 28;
	poly_space.resize( poly_size );
    
	double x = center[0];
	double y = center[1];

	// 1st order.
	poly_space[0] = 1.0;
	poly_space[1] = x;
	poly_space[2] = y;

	// 2nd order.
	double x2 = x*x;
	double y2 = y*y;
	poly_space[3] = x2;
	poly_space[4] = x*y;
	poly_space[5] = y2;

	// 3rd order.
	double x3 = x2*x;
	double y3 = y2*y;
	poly_space[6] = x3;
	poly_space[7] = x2*y;
	poly_space[8] = x*y2;
	poly_space[9] = y3;

	// 4th order.
	double x4 = x3*x;
	double y4 = y3*y;
	poly_space[10] = x4;
	poly_space[11] = x3*y;
	poly_space[12] = x2*y2;
	poly_space[13] = x*y3;
	poly_space[14] = y4;

	// 5th order.
	double x5 = x4*x;
	double y5 = y4*y;
	poly_space[15] = x5;
	poly_space[16] = x4*y;
	poly_space[17] = x3*y2;
	poly_space[18] = x3*y2;
	poly_space[19] = x*y4;
	poly_space[20] = y5;

	// 6th order.
	double x6 = x5*x;
	double y6 = y5*y;
	poly_space[21] = x6;
	poly_space[22] = x5*y;
	poly_space[23] = x4*y2;
	poly_space[24] = x3*y3;
	poly_space[25] = x2*y4;
	poly_space[26] = x*y5;
	poly_space[27] = y6;
    }

    // 3 dimensions.
    else if ( 3 == DIM )
    {
	int poly_size = 56;
	poly_space.resize( poly_size );
    
	double x = center[0];
	double y = center[1];
	double z = center[2];

	// 1st order.
	poly_space[0] = 1.0;
	poly_space[1] = x;
	poly_space[2] = y;
	poly_space[3] = z;

	// 2nd order.
	double x2 = x*x;
	double y2 = y*y;
	double z2 = z*z;
	poly_space[4] = x*y;
	poly_space[5] = x*z;
	poly_space[6] = y*z;
	poly_space[7] = x2;
	poly_space[8] = y2;
	poly_space[9] = z2;

	// 3rd order.
	double x3 = x2*x;
	double y3 = y2*y;
	double z3 = z2*z;
	poly_space[10] = x*y*z;
	poly_space[11] = x2*y;
	poly_space[12] = x2*z;
	poly_space[13] = x*y2;
	poly_space[14] = y2*z;
	poly_space[15] = x*z2;
	poly_space[16] = y*z2;
	poly_space[17] = x3;
	poly_space[18] = y3;
	poly_space[19] = z3;

	// 4th order.
	double x4 = x3*x;
	double y4 = y3*y;
	double z4 = z3*z;
	poly_space[20] = x2*y*z;
	poly_space[21] = x*y2*z;
	poly_space[22] = x*y*z2;
	poly_space[23] = x3*y;
	poly_space[24] = x3*z;
	poly_space[25] = x*y3;
	poly_space[26] = y3*z;
	poly_space[27] = x*z3;
	poly_space[28] = y*z3;
	poly_space[29] = x2*y2;
	poly_space[30] = x2*z2;
	poly_space[31] = y2*z2;
	poly_space[32] = x4;
	poly_space[33] = y4;
	poly_space[34] = z4;

	// 5th order.
	double x5 = x4*x;
	double y5 = y4*y;
	double z5 = z4*z;
	poly_space[35] = x3*y*z;
	poly_space[36] = x*y3*z;
	poly_space[37] = x*y*z3;
	poly_space[38] = x2*y2*z;
	poly_space[39] = x2*y*z2;
	poly_space[40] = x*y2*z2;
	poly_space[41] = x3*y2;
	poly_space[42] = x3*z2;
	poly_space[43] = x2*y3;
	poly_space[44] = y3*z2;
	poly_space[45] = x2*z3;
	poly_space[46] = y2*z3;
	poly_space[47] = x4*y;
	poly_space[48] = x4*z;
	poly_space[49] = x*y4;
	poly_space[50] = y4*z;
	poly_space[51] = x*z4;
	poly_space[52] = y*z4;
	poly_space[53] = x5;
	poly_space[54] = y5;
	poly_space[55] = z5;
    }

    return poly_space;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given an SPD matrix, compute the LDLT factorization.
 */
template<class Basis, class GO, int DIM>
void OrthogonalPolynomialMLSProblem<Basis,GO,DIM>::factorLDLT( 
    Teuchos::SerialDenseMatrix<int,double>& A ) const
{
    DTK_CHECK( A.numRows() == A.numCols() );

    int n = A.numRows();
    Teuchos::Array<double> v(n);
    for ( int j = 0; j < n; ++j )
    {
	for ( int i = 0; i < j; ++i )
	{
	    v[i] = A(j,i)*A(i,i);
	}
	v[j] = A(j,j);
	for ( int i = 0; i < j; ++i )
	{
	    v[j] -= A(j,i)*v[i];
	}
	A(j,j) = v[j];
	for ( int i = j+1; i < n; ++i )
	{
	    for ( int k = 0; k < j; ++k )
	    {
		A(i,j) -= A(i,k)*v[k];
	    }
	    DTK_CHECK( v[j] != 0.0 );
	    A(i,j) /= v[j];
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_ORTHOGONALPOLYNOMIALMLSPROBLEM_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_OrthogonalPolynomialMLSProblem_impl.hpp
//---------------------------------------------------------------------------//

