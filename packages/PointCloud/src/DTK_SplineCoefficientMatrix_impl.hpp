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
template<class Basis, class GO, int DIM>
SplineCoefficientMatrix<Basis,GO,DIM>::SplineCoefficientMatrix(
    Teuchos::RCP<const Tpetra::Map<int,GO> >& operator_map,
    const Teuchos::ArrayView<const double>& source_centers,
    const Teuchos::ArrayView<const GO>& source_center_gids,
    const Teuchos::ArrayView<const double>& dist_source_centers,
    const Teuchos::ArrayView<const GO>& dist_source_center_gids,
    const Teuchos::RCP<SplineInterpolationPairing<DIM> >& source_pairings,
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

    // Create the P^T matrix.
    int offset = DIM + 1;
    unsigned lid_offset = operator_map->getComm()->getRank() ? 0 : offset;
    Teuchos::RCP<Tpetra::MultiVector<double,int,GO> > P_trans_vec = 
	Tpetra::createMultiVector<double,int,GO>( operator_map, offset );
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > P_trans_view = 
	P_trans_vec->get2dViewNonConst();
    int di = 0; 
    for ( unsigned i = lid_offset; i < num_source_centers + lid_offset; ++i )
    {
	P_trans_view[0][i] = 1.0;
	di = DIM*(i-lid_offset);
	for ( int d = 0; d < DIM; ++d )
	{
	    P_trans_view[d+1][i] = source_centers[di+d];
	}
    }
    d_P_trans =Teuchos::rcp( new PolynomialMatrix<GO>(P_trans_vec) );

    // Create the M matrix.
    Teuchos::ArrayRCP<std::size_t> entries_per_row;
    if ( 0 == operator_map->getComm()->getRank() )
    {
	entries_per_row = Teuchos::ArrayRCP<std::size_t>(
	    num_source_centers + offset, 0 );
	Teuchos::ArrayRCP<std::size_t> children_per_parent =
	    source_pairings->childrenPerParent();
	std::copy( children_per_parent.begin(), children_per_parent.end(),
		   entries_per_row.begin() + offset );
    }
    else
    {
	entries_per_row = source_pairings->childrenPerParent();
    }
    d_M = Teuchos::rcp( new Tpetra::CrsMatrix<double,int,GO>(
			    operator_map,
			    entries_per_row, Tpetra::StaticProfile) );
    Teuchos::Array<GO> M_indices;
    Teuchos::Array<double> values;
    int dj = 0;
    Teuchos::ArrayView<const unsigned> source_neighbors;
    double dist = 0.0;
    for ( unsigned i = 0; i < num_source_centers; ++i )
    {
    	di = DIM*i;
	source_neighbors = source_pairings->childCenterIds( i );
	M_indices.resize( source_neighbors.size() );
	values.resize( source_neighbors.size() );
    	for ( unsigned j = 0; j < source_neighbors.size(); ++j )
    	{
	    dj = DIM*source_neighbors[j];
	    M_indices[j] = 
		dist_source_center_gids[ source_neighbors[j] ];

	    dist = EuclideanDistance<DIM>::distance(
		&source_centers[di], &dist_source_centers[dj] );

    	    values[j] = BP::evaluateValue( basis, dist );
    	}

	d_M->insertGlobalValues( source_center_gids[i], M_indices(), values() );
    }
    d_M->fillComplete();

    DTK_ENSURE( d_M->isFillComplete() );
}

//---------------------------------------------------------------------------//
// Apply operation. The operator is symmetric and therefore the transpose
// apply is the same as the forward apply.
template<class Basis, class GO, int DIM>
void SplineCoefficientMatrix<Basis,GO,DIM>::apply(
    const Tpetra::MultiVector<double,int,GO> &X,
    Tpetra::MultiVector<double,int,GO> &Y,
    Teuchos::ETransp mode,
    double alpha,
    double beta ) const
{
    // Make a work vector.
    Tpetra::MultiVector<double,int,GO> work( X );

    // Apply P^T.
    d_P_trans->apply( X, Y, Teuchos::NO_TRANS, alpha, beta );

    // Apply P.
    d_P_trans->apply( X, work, Teuchos::TRANS, alpha, beta );
    Y.update( 1.0, work, 1.0 );

    // Apply M.
    d_M->apply( X, work, Teuchos::NO_TRANS, alpha, beta );
    Y.update( 1.0, work, 1.0 );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_SPLINECOEFFICIENTMATRIX_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineCoefficientMatrix_impl.hpp
//---------------------------------------------------------------------------//
