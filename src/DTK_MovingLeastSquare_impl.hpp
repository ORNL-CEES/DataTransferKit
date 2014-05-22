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
 * \file   DTK_MovingLeastSquare_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Moving least square interpolator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOVINGLEASTSQUARE_IMPL_HPP
#define DTK_MOVINGLEASTSQUARE_IMPL_HPP

#include "DTK_DBC.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_MultiVector.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param comm The parallel communicator over which the perform the
 * interpolation.
 */
//---------------------------------------------------------------------------//
template<class Basis, class GO, int DIM>
MovingLeastSquare<Basis,GO,DIM>::MovingLeastSquare( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm )
    : d_comm( comm )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Set the interpolation problem.
 *
 * \param source_centers The cartesian coordinates of the source
 * centers. Coordinates must be interleaved
 * (x0,y0,z0,x1,y1,z1,...,xN,yN,zN). If there are no source centers on this
 * process then provide an array view of size 0.
 *
 * \param target_centers The cartesian coordinates of the target
 * centers. Coordinates must be interleaved
 * (x0,y0,z0,x1,y1,z1,...,xN,yN,zN). If there are no target centers on this
 * process then provide an array view of size 0.
 *
 * \brief radius The radius over which the basis functions will be
 * defined. Must be greater than 0.
 */
template<class Basis, class GO, int DIM>
void MovingLeastSquare<Basis,GO,DIM>::setProblem(
    const Teuchos::ArrayView<const double>& source_centers,
    const Teuchos::ArrayView<const double>& target_centers,
    const double radius )
{
    DTK_REQUIRE( 0 == source_centers.size() % DIM );
    DTK_REQUIRE( 0 == target_centers.size() % DIM );
    DTK_REQUIRE( radius > 0.0 );

    // Build the interpolation matrix.
    buildInterpolationMatrix( source_centers, target_centers, radius );

    DTK_ENSURE( Teuchos::nonnull(d_H) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given a set of scalar values at the given source centers in the
 *  source decomposition, interpolate them onto the target centers in the
 *  target decomposition.
 *
 * \param source_data View of the source data defined at the source
 * centers. The data must be blocked by dimension. If there is no data on this
 * process then the view must be of size 0.
 *
 * \param num_source_dims Number of source data dimensions. Must be the same
 * as the number of target data dimensions.
 *
 * \param source_lda The stride of the source vectors. Must be equal to the
 * number of source centers.
 *
 * \param target_data View of the target data defined at the target
 * centers. The data must be blocked by dimension. If there is no data on this
 * process then the view must be of size 0.
 *
 * \param num_target_dims Number of target data dimensions. Must be the same
 * as the number of target data dimensions.
 *
 * \param target_lda The stride of the target vectors. Must be equal to the
 * number of target centers.
 *
 * \param max_solve_iterations Maximum number of linear solver iterations
 * allowed in the interpolation solution.
 *
 * \param solve_convergence_tolerance Linear solver convergence tolerance for
 * the interpolation solution.
 */
template<class Basis, class GO, int DIM>
void MovingLeastSquare<Basis,GO,DIM>::interpolate( 
    const Teuchos::ArrayView<const double>& source_data,
    const int num_source_dims,
    const int source_lda,
    const Teuchos::ArrayView<double>& target_data,
    const int num_target_dims,
    const int target_lda ) const
{
    DTK_REQUIRE( num_source_dims == num_target_dims );
    DTK_REQUIRE( source_data.size() == source_lda * num_source_dims );
    DTK_REQUIRE( target_data.size() == target_lda * num_target_dims );

    // Build the source vector.
    Teuchos::ArrayRCP<double> source_data_view(
	const_cast<double*>(source_data.getRawPtr()), 0, 
	source_data.size(), false );
    Teuchos::RCP<Tpetra::MultiVector<double,int,GO> >
	source_vec = Tpetra::createMultiVectorFromView( 
	    d_source_map, source_data_view, source_lda, num_source_dims );

    // Build the target vector.
    Teuchos::ArrayRCP<double> target_data_view =
	Teuchos::arcpFromArrayView(target_data);
    Teuchos::RCP<Tpetra::MultiVector<double,int,GO> >
	target_vec = Tpetra::createMultiVectorFromView( 
	    d_target_map, target_data_view, target_lda, num_target_dims );

    // Apply the interpolation matrix.
    d_H->apply( *source_vec, *target_vec );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the interpolation matrix.
 */
template<class Basis, class GO, int DIM>
void MovingLeastSquare<Basis,GO,DIM>::buildInterpolationMatrix(
    const Teuchos::ArrayView<const double>& source_centers,
    const Teuchos::ArrayView<const double>& target_centers,
    const double radius )
{
    // Build the source map.
    GO local_num_src = source_centers.size() / DIM;
    GO global_num_src = 0;
    Teuchos::reduceAll<int,GO>( *d_comm, 
				Teuchos::REDUCE_SUM,
				local_num_src,
				Teuchos::Ptr<GO>(&global_num_src) );
    d_source_map = Tpetra::createContigMap<int,GO>( global_num_src,
						    local_num_src,
						    d_comm );
    Teuchos::ArrayView<const GO> source_gids = 
	d_source_map->getNodeElementList();

    // Build the target map.
    GO local_num_tgt = target_centers.size() / DIM;
    GO global_num_tgt = 0;
    Teuchos::reduceAll<int,GO>( *d_comm, 
				Teuchos::REDUCE_SUM,
				local_num_tgt,
				Teuchos::Ptr<GO>(&global_num_tgt) );
    d_target_map = Tpetra::createContigMap<int,GO>( global_num_tgt,
						    local_num_tgt,
						    d_comm );
    Teuchos::ArrayView<const GO> target_gids = 
	d_target_map->getNodeElementList();

    // Build the basis.
    Teuchos::RCP<Basis> basis = BP::create( radius );

    // Gather the source centers that are within a radius of the target
    // centers on this proc.
    Teuchos::Array<double> dist_sources;
    CenterDistributor<DIM> distributor( 
	d_comm, source_centers, target_centers, radius, dist_sources );

    // Gather the global ids of the source centers that are within a radius of
    // the target centers on this proc.
    Teuchos::Array<GO> dist_source_gids( distributor.getNumImports() );
    distributor.distribute( source_gids, dist_source_gids() );

    // Build the source/target pairings.
    SplineInterpolationPairing<DIM> pairings( 
	dist_sources, target_centers, radius );

    // Build the interpolation matrix.
    d_H = Tpetra::createCrsMatrix<double,int,GO>( d_target_map );
    Teuchos::ArrayView<const double> target_view;
    Teuchos::Array<GO> indices;
    Teuchos::ArrayView<const double> values;
    Teuchos::ArrayView<const unsigned> pair_gids;
    for ( int i = 0; i < local_num_tgt; ++i )
    {
	// If there is no support for this target center then do not build a
	// local basis.
	if ( 0 < pairings.childCenterIds(i).size() )
	{
	    // Get a view of this target center.
	    target_view = target_centers(i*DIM,DIM);

	    // Build the local interpolation problem. 
	    LocalMLSProblem<Basis,GO,DIM> local_problem(
		target_view, pairings.childCenterIds(i),
		dist_sources, *basis );

	    // Populate the interpolation matrix row.
	    values = local_problem.shapeFunction();
	    indices.resize( values.size() );
	    pair_gids = pairings.childCenterIds(i);
	    for ( unsigned j = 0; j < values.size(); ++j )
	    {
		indices[j] = dist_source_gids[ pair_gids[j] ];
	    }
	    d_H->insertGlobalValues( target_gids[i], indices(), values );
	}
    }

    d_H->fillComplete( d_source_map, d_target_map );
    DTK_ENSURE( d_H->isFillComplete() );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MOVINGLEASTSQUARE_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_MovingLeastSquare_impl.hpp
//---------------------------------------------------------------------------//

