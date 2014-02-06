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
 * \file   DTK_SplineOperatorA_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Spline transformation operator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SPLINEOPERATORA_IMPL_HPP
#define DTK_SPLINEOPERATORA_IMPL_HPP

#include "DTK_DBC.hpp"
#include "DTK_EuclideanDistance.hpp"

#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Build the transformation matrix.
 */
template<class Basis, class GO, int DIM>
Teuchos::RCP<Tpetra::CrsMatrix<double,int,GO> > 
SplineOperatorA<Basis,GO,DIM>::create(
    Teuchos::RCP<const Tpetra::Map<int,GO> >& domain_map,
    Teuchos::RCP<const Tpetra::Map<int,GO> >& range_map,
    const Teuchos::ArrayView<const double>& target_centers,
    const Teuchos::ArrayView<const GO>& target_center_gids,
    const Teuchos::ArrayView<const double>& dist_source_centers,
    const Teuchos::ArrayView<const GO>& dist_source_center_gids,
    const Teuchos::RCP<SplineInterpolationPairing<DIM> >& target_pairings,
    const Basis& basis,
    const double alpha )
{
    DTK_CHECK( 0 == target_centers.size() % DIM );
    DTK_CHECK( target_centers.size() / DIM == 
		    target_center_gids.size() );
    DTK_CHECK( 0 == dist_source_centers.size() % DIM );
    DTK_CHECK( dist_source_centers.size() / DIM == 
		    dist_source_center_gids.size() );

    // Create the matrix.
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,GO> > A =
	Tpetra::createCrsMatrix<double,int,GO>( domain_map );

    unsigned num_target_centers = target_center_gids.size();

    int offset = DIM+1;
    int di = 0;
    int dj = 0;
    int jpoffset = 0;
    Teuchos::ArrayView<const unsigned> target_neighbors;
    Teuchos::Array<double> values;
    Teuchos::Array<GO> indices;
    double dist = 0.0;
    for ( unsigned i = 0; i < num_target_centers; ++i )
    {
	di = DIM*i;

	// Get the source points neighboring this target point.
	target_neighbors = target_pairings->childCenterIds( i );
	values.resize( offset + target_neighbors.size() );
	indices.resize( offset + target_neighbors.size() );

	// 1's column.
	indices[0] = 0;
	values[0] = 1.0;

	// Add the coordinates.
	for ( int d = 0; d < DIM; ++d )
	{
	    indices[d+1] = d+1;
	    values[d+1] = target_centers[di+d];
	}

	// Add the local basis contributions.
    	for ( unsigned j = 0; j < target_neighbors.size(); ++j )
    	{
	    dj = DIM*target_neighbors[j];
	    jpoffset = j+offset;

	    indices[jpoffset] = 
		dist_source_center_gids[ target_neighbors[j] ];

	    dist = EuclideanDistance<DIM>::distance(
		&target_centers[di], &dist_source_centers[dj] );

    	    values[jpoffset] = BP::evaluateValue( basis, dist );
    	    if ( alpha > 0.0 )
    	    {
    		values[jpoffset] += alpha * BP::evaluateGradient( basis, dist );
    	    }
    	}

	A->insertGlobalValues( target_center_gids[i], indices(), values() );
    }

    A->fillComplete( range_map, domain_map );
    DTK_ENSURE( A->isFillComplete() );

    return A;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_SPLINEOPERATORA_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineOperatorA_impl.hpp
//---------------------------------------------------------------------------//

