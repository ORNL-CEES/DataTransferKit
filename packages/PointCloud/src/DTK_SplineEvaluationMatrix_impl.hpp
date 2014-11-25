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
 * \file   DTK_SplineEvaluationMatrix_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Spline transformation operator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SPLINEEVALUATIONMATRIX_IMPL_HPP
#define DTK_SPLINEEVALUATIONMATRIX_IMPL_HPP

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
SplineEvaluationMatrix<Basis,DIM>::SplineEvaluationMatrix(
    const Teuchos::RCP<const Tpetra::Map<int,std::size_t> >& domain_map,
    const Teuchos::RCP<const Tpetra::Map<int,std::size_t> >& range_map,
    const Teuchos::ArrayView<const double>& target_centers,
    const Teuchos::ArrayView<const std::size_t>& target_center_gids,
    const Teuchos::ArrayView<const double>& dist_source_centers,
    const Teuchos::ArrayView<const std::size_t>& dist_source_center_gids,
    const SplineInterpolationPairing<DIM>& target_pairings,
    const Basis& basis )
{
    DTK_CHECK( 0 == target_centers.size() % DIM );
    DTK_CHECK( target_centers.size() / DIM == 
		    target_center_gids.size() );
    DTK_CHECK( 0 == dist_source_centers.size() % DIM );
    DTK_CHECK( dist_source_centers.size() / DIM == 
		    dist_source_center_gids.size() );

    // Get the number of target centers.
    unsigned num_target_centers = target_center_gids.size();

    // Create the Q matrix.
    int offset = DIM + 1;
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > Q_vec = 
	Tpetra::createMultiVector<double,int,std::size_t>( range_map, offset );
    int di = 0; 
    for ( unsigned i = 0; i < num_target_centers; ++i )
    {
	Q_vec->replaceGlobalValue( target_center_gids[i], 0, 1.0 );
	di = DIM*i;
	for ( int d = 0; d < DIM; ++d )
	{
	    Q_vec->replaceGlobalValue( 
		target_center_gids[i], d+1, target_centers[di+d] );
	}
    }
    d_Q = Teuchos::rcp( 
	new PolynomialMatrix<std::size_t>(Q_vec,domain_map,range_map) );

    // Create the N matrix.
    Teuchos::ArrayRCP<std::size_t> children_per_parent =
	target_pairings.childrenPerParent();
    std::size_t max_entries_per_row = *std::max_element( 
	children_per_parent.begin(), children_per_parent.end() );
    d_N = Teuchos::rcp( new Tpetra::CrsMatrix<double,int,std::size_t>( 
			    range_map, max_entries_per_row) );
    Teuchos::Array<std::size_t> N_indices( max_entries_per_row );
    Teuchos::Array<double> values( max_entries_per_row );
    int dj = 0;
    Teuchos::ArrayView<const unsigned> target_neighbors;
    double dist = 0.0;
    int ntn = 0;
    for ( unsigned i = 0; i < num_target_centers; ++i )
    {
	di = DIM*i;

	// Get the source points neighboring this target point.
	target_neighbors = target_pairings.childCenterIds( i );
	ntn = target_neighbors.size();

	// Add the local basis contributions.
    	for ( unsigned j = 0; j < ntn; ++j )
    	{
	    dj = DIM*target_neighbors[j];

	    N_indices[j] = 
		dist_source_center_gids[ target_neighbors[j] ];

	    dist = EuclideanDistance<DIM>::distance(
		&target_centers[di], &dist_source_centers[dj] );

    	    values[j] = BP::evaluateValue( basis, dist );
    	}

	d_N->insertGlobalValues( 
	    target_center_gids[i], N_indices(0,ntn), values(0,ntn) );
    }
    d_N->fillComplete( domain_map, range_map );

    DTK_ENSURE( d_N->isFillComplete() );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_SPLINEEVALUATIONMATRIX_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineEvaluationMatrix_impl.hpp
//---------------------------------------------------------------------------//

