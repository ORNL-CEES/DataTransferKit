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
 * \file   DTK_SplineCoefficientMatrix_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Spline coefficient matrix.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SPLINECOEFFICIENTMATRIX_IMPL_HPP
#define DTK_SPLINECOEFFICIENTMATRIX_IMPL_HPP

#include "DTK_DBC.hpp"
#include "DTK_EuclideanDistance.hpp"

#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Basis,int DIM>
SplineCoefficientMatrix<Basis,DIM>::SplineCoefficientMatrix(
    const Teuchos::RCP<const Tpetra::Map<int,std::size_t> >& operator_map,
    const Teuchos::ArrayView<const double>& source_centers,
    const Teuchos::ArrayView<const std::size_t>& source_center_gids,
    const Teuchos::ArrayView<const double>& dist_source_centers,
    const Teuchos::ArrayView<const std::size_t>& dist_source_center_gids,
    const SplineInterpolationPairing<DIM>& source_pairings,
    const Basis& basis )
{
    DTK_CHECK( 0 == source_centers.size() % DIM );
    DTK_CHECK( source_centers.size() / DIM == 
		    source_center_gids.size() );
    DTK_CHECK( 0 == dist_source_centers.size() % DIM );
    DTK_CHECK( dist_source_centers.size() / DIM == 
		    dist_source_center_gids.size() );

    // Get the number of source centers.
    unsigned num_source_centers = source_center_gids.size();

    // Create the P matrix.
    int offset = DIM + 1;
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > P_vec = 
	Tpetra::createMultiVector<double,int,std::size_t>( operator_map, offset );
    int di = 0; 
    for ( unsigned i = 0; i < num_source_centers; ++i )
    {
	P_vec->replaceGlobalValue( source_center_gids[i], 0, 1.0 );
	di = DIM*i;
	for ( int d = 0; d < DIM; ++d )
	{
	    P_vec->replaceGlobalValue( 
		source_center_gids[i], d+1, source_centers[di+d] );
	}
    }
    d_P =Teuchos::rcp(
	new PolynomialMatrix<std::size_t>(P_vec,operator_map,operator_map) );

    // Create the M matrix.
    Teuchos::ArrayRCP<std::size_t> children_per_parent =
	source_pairings.childrenPerParent();
    std::size_t max_entries_per_row = *std::max_element( 
	children_per_parent.begin(), children_per_parent.end() );
    d_M = Teuchos::rcp( new Tpetra::CrsMatrix<double,int,std::size_t>(
			    operator_map, max_entries_per_row) );
    Teuchos::Array<std::size_t> M_indices( max_entries_per_row );
    Teuchos::Array<double> values( max_entries_per_row );
    int dj = 0;
    Teuchos::ArrayView<const unsigned> source_neighbors;
    double dist = 0.0;
    int nsn = 0;
    for ( unsigned i = 0; i < num_source_centers; ++i )
    {
	// Get the source points neighboring this source point.
    	di = DIM*i;
	source_neighbors = source_pairings.childCenterIds( i );
	nsn = source_neighbors.size();

	// Add the local basis contributions.
    	for ( int j = 0; j < nsn; ++j )
    	{
	    dj = DIM*source_neighbors[j];
	    M_indices[j] = 
		dist_source_center_gids[ source_neighbors[j] ];

	    dist = EuclideanDistance<DIM>::distance(
		&source_centers[di], &dist_source_centers[dj] );

    	    values[j] = BP::evaluateValue( basis, dist );
    	}
	d_M->insertGlobalValues( 
	    source_center_gids[i], M_indices(0,nsn), values(0,nsn) );
    }
    d_M->fillComplete();

    DTK_ENSURE( d_M->isFillComplete() );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_SPLINECOEFFICIENTMATRIX_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineCoefficientMatrix_impl.hpp
//---------------------------------------------------------------------------//
