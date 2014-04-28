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
 * \file   DTK_SplineOperatorC_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Spline interpolation operator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SPLINEOPERATORC_IMPL_HPP
#define DTK_SPLINEOPERATORC_IMPL_HPP

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
SplineOperatorC<Basis,GO,DIM>::SplineOperatorC(
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
    d_P_trans = Teuchos::rcp( 
	new Tpetra::CrsMatrix<double,int,GO>( 
	    operator_map, 1 + DIM, Tpetra::StaticProfile) );

    int offset = DIM+1;
    int di = 0; 
    Teuchos::Array<GO> indices(offset,0);
    Teuchos::Array<double> values(offset,1);
    for ( int i = 1; i < offset; ++i )
    {
	indices[i] = i;
    }
    for ( unsigned i = 0; i < num_source_centers; ++i )
    {
	di = DIM*i;

	for ( int d = 0; d < DIM; ++d )
	{
	    values[d+1] = source_centers[di+d];
	}

	d_P_trans->insertGlobalValues( 
	    source_center_gids[i], indices(), values() );
    }
    d_P_trans->fillComplete();

    // Create the M matrix.
    d_M = Tpetra::createCrsMatrix<double,int,GO>( operator_map );
    int dj = 0;
    Teuchos::ArrayView<const unsigned> source_neighbors;
    double dist = 0.0;
    for ( unsigned i = 0; i < num_source_centers; ++i )
    {
    	di = DIM*i;
	source_neighbors = source_pairings->childCenterIds( i );
	indices.resize( source_neighbors.size() );
	values.resize( source_neighbors.size() );
    	for ( unsigned j = 0; j < source_neighbors.size(); ++j )
    	{
	    dj = DIM*source_neighbors[j];
	    indices[j] = 
		dist_source_center_gids[ source_neighbors[j] ];

	    dist = EuclideanDistance<DIM>::distance(
		&source_centers[di], &dist_source_centers[dj] );

    	    values[j] = BP::evaluateValue( basis, dist );
    	}

	d_M->insertGlobalValues( source_center_gids[i], indices(), values() );
    }
    d_M->fillComplete();

    DTK_ENSURE( d_P_trans->isFillComplete() );
    DTK_ENSURE( d_M->isFillComplete() );
}

//---------------------------------------------------------------------------//
// Apply operation. The operator is symmetric and therefore the transpose
// apply is the same as the forward apply.
template<class Basis, class GO, int DIM>
void SplineOperatorC<Basis,GO,DIM>::apply(
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

#endif // end DTK_SPLINEOPERATORC_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineOperatorC_impl.hpp
//---------------------------------------------------------------------------//

