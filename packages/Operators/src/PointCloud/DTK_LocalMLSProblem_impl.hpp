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
 * \file   DTK_LocalMLSProblem_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Local moving least square problem.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_LOCALMLSPROBLEM_IMPL_HPP
#define DTK_LOCALMLSPROBLEM_IMPL_HPP

#include <limits>
#include <functional>
#include <algorithm>

#include "DTK_LocalMLSProblem.hpp"
#include "DTK_RadialBasisPolicy.hpp"
#include "DTK_EuclideanDistance.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_LAPACK.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Basis,int DIM>
LocalMLSProblem<Basis,DIM>::LocalMLSProblem(
    const Teuchos::ArrayView<const double>& target_center,
    const Teuchos::ArrayView<const unsigned>& source_lids,
    const Teuchos::ArrayView<const double>& source_centers,
    const Basis& basis,
    const double radius )
    : d_shape_function( source_lids.size() )
{
    DTK_REQUIRE( 0 == source_centers.size() % DIM );
    DTK_REQUIRE( 0 == target_center.size() % DIM );

    // Number of source centers supporting this target center.
    int num_sources = source_lids.size();
    int poly_size = 0;
    Teuchos::SerialDenseMatrix<int,double> P;
    Teuchos::SerialDenseVector<int,double> target_poly;

    // Make Phi.
    Teuchos::SerialDenseMatrix<int,double> phi( num_sources, num_sources );
    Teuchos::ArrayView<const double> source_center_view;
    double dist = 0.0;
    for ( int i = 0; i < num_sources; ++i )
    {
        source_center_view = source_centers(DIM*source_lids[i],DIM);
        dist = EuclideanDistance<DIM>::distance(
            target_center.getRawPtr(), source_center_view.getRawPtr() );
        phi(i,i) = BP::evaluateValue( basis, radius, dist );
    }

    // Make P.
    Teuchos::Array<int> poly_ids(1,0);
    int poly_id = 0;
    P.reshape( num_sources, poly_id+1 );
    for ( int i = 0; i < num_sources; ++i )
    {
        source_center_view = source_centers(DIM*source_lids[i],DIM);
        P(i,poly_id) = polynomialCoefficient( 0, source_center_view );
    }
    ++poly_id;

    // Add polynomial columns until we are full rank.
    bool full_rank = false;
    int num_poly = 10;
    int total_poly = std::min( num_poly, num_sources );
    for ( int j = 1; j < total_poly; ++j )
    {
        // Add the next column.
        P.reshape( num_sources, poly_id+1 );
        for ( int i = 0; i < num_sources; ++i )
        {
            source_center_view = source_centers(DIM*source_lids[i],DIM);
            P(i,poly_id) = polynomialCoefficient( j, source_center_view );
        }

        // Check for rank deficiency.
        full_rank = isFullRank( P );

        // If we are full rank, add this coefficient.
        if ( full_rank )
        {
            poly_ids.push_back( j );
            ++poly_id;
        }

        // If we are rank deficient, remove the last column.
        else
        {
            P.reshape( num_sources, poly_id );
        }
    }

    // Make p.
    poly_size = poly_ids.size();
    target_poly.resize( poly_size );
    for ( int i = 0; i < poly_size; ++i )
    {
        target_poly(i) = polynomialCoefficient( poly_ids[i], target_center );
    }

    // Construct b.
    Teuchos::SerialDenseMatrix<int,double> b( poly_size, num_sources );
    b.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, P, phi, 0.0 );

    // Construct the A matrix.
    Teuchos::SerialDenseMatrix<int,double> A( poly_size, poly_size );
    {
        // Build A.
        Teuchos::SerialDenseMatrix<int,double> work( num_sources, poly_size );
        work.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, phi, P, 0.0 );
        A.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, P, work, 0.0 );
    }

    // Apply the inverse of the A matrix to b.
    Teuchos::LAPACK<int,double> lapack;
    double A_rcond = std::numeric_limits<double>::epsilon();
    Teuchos::Array<double> work( 4 * A.numCols() );
    Teuchos::SerialDenseVector<int,double> s( poly_size );
    int rank = 0;
    int info = 0;

    // Estimate the reciprocal of the condition number.
    Teuchos::Array<int> ipiv( std::min(A.numRows(),A.numCols()) );
    Teuchos::SerialDenseMatrix<int,double> LU_A( A );
    lapack.GETRF( LU_A.numRows(), LU_A.numCols(), LU_A.values(),
                  LU_A.numRows(), ipiv.getRawPtr(), &info );
    DTK_CHECK( 0 == info );

    Teuchos::Array<int> iwork( A.numCols() );
    lapack.GECON( '1', LU_A.numCols(), LU_A.values(),
                  LU_A.numRows(), A.normOne(), &A_rcond,
                  work.getRawPtr(), iwork.getRawPtr(), &info );
    DTK_CHECK( 0 == info );

    // Get the optimal work size.
    lapack.GELSS( A.numRows(), A.numCols(), b.numCols(),
                  A.values(), A.numRows(),
                  b.values(), b.numRows(), s.values(),
                  A_rcond, &rank, work.getRawPtr(), -1, &info );
    DTK_CHECK( 0 == info );

    // Apply the inverse of A to b.
    work.resize( work[0] );
    lapack.GELSS( A.numRows(), A.numCols(), b.numCols(),
                  A.values(), A.numRows(),
                  b.values(), b.numRows(), s.values(),
                  A_rcond, &rank, work.getRawPtr(), work.size(), &info );
    DTK_CHECK( 0 == info );

    // Construct the basis.
    Teuchos::SerialDenseMatrix<int,double> shape_matrix(
        Teuchos::View, d_shape_function.getRawPtr(),
        1, 1, d_shape_function.size() );
    shape_matrix.multiply( Teuchos::TRANS, Teuchos::NO_TRANS,
                           1.0, target_poly, b, 0.0 );
}

//---------------------------------------------------------------------------//
// Get a polynomial coefficient.
template<class Basis,int DIM>
double LocalMLSProblem<Basis,DIM>::polynomialCoefficient(
    const int coeff, const Teuchos::ArrayView<const double>& center ) const
{
    switch( coeff )
    {
        // Linear.
        case( 0 ):
            return 1.0;
            break;
        case( 1 ):
            return center[0];
            break;
        case( 2 ):
            return center[1];
            break;
        case( 3 ):
            return center[2];
            break;

        // Quadratic
        case( 4 ):
            return center[0]*center[1];
            break;
        case( 5 ):
            return center[0]*center[2];
            break;
        case( 6 ):
            return center[1]*center[2];
            break;
        case( 7 ):
            return center[0]*center[0];
            break;
        case( 8 ):
            return center[1]*center[1];
            break;
        case( 9 ):
            return center[2]*center[2];
            break;
    }
    return 0.0;
}

//---------------------------------------------------------------------------//
// Check if a matrix is full rank.
template<class Basis,int DIM>
bool LocalMLSProblem<Basis,DIM>::isFullRank(
    const Teuchos::SerialDenseMatrix<int,double>& matrix ) const
{
    // Copy the matrix.
    Teuchos::SerialDenseMatrix<int,double> A = matrix;

    // Determine the full rank.
    int full_rank = std::min( A.numRows(), A.numCols() );

    // Compute the singular value decomposition.
    Teuchos::LAPACK<int,double> lapack;
    Teuchos::Array<double> S( full_rank );
    Teuchos::SerialDenseMatrix<int,double> U( A.numRows(), A.numRows() );
    Teuchos::SerialDenseMatrix<int,double> VT( A.numCols(), A.numCols() );
    Teuchos::Array<double> work( full_rank );
    Teuchos::Array<double> rwork( full_rank );
    int info = 0;
    lapack.GESVD( 'A', 'A',
                  A.numRows(), A.numCols(), A.values(), A.numRows(),
                  S.getRawPtr(),
                  U.values(), U.numRows(),
                  VT.values(), VT.numRows(),
                  work.getRawPtr(), -1, rwork.getRawPtr(), &info );
    DTK_CHECK( 0 == info );

    work.resize( work[0] );
    rwork.resize( work.size() );
    lapack.GESVD( 'A', 'A',
                  A.numRows(), A.numCols(), A.values(), A.numRows(),
                  S.getRawPtr(),
                  U.values(), U.numRows(),
                  VT.values(), VT.numRows(),
                  work.getRawPtr(), work.size(), rwork.getRawPtr(), &info );
    DTK_CHECK( 0 == info );

    // Check the singular values. If they are greater than delta they count.
    double epsilon = std::numeric_limits<double>::epsilon();
    double delta = S[0]*epsilon;
    int rank = std::count_if( S.begin(), S.end(),
                              [=](double s){ return (s > delta); } );

    return ( rank == full_rank );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_LOCALMLSPROBLEM_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_LocalMLSProblem_impl.hpp
//---------------------------------------------------------------------------//

