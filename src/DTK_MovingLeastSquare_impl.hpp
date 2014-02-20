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
#include "DTK_CenterDistributor.hpp"
#include "DTK_SplineInterpolationPairing.hpp"
#include "DTK_SplineOperatorC.hpp"
#include "DTK_SplineOperatorA.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Map.hpp>

#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosLinearProblem.hpp>

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
MovingLeastSquares<Basis,GO,DIM>::MovingLeastSquare( 
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
 *
 * \brief alpha Derivative contribution parameter. Must be greater than or
 * equal to zero. A value of zero will give no derivative contributions to the
 * interpolation.
 */
template<class Basis, class GO, int DIM>
void MovingLeastSquare<Basis,GO,DIM>::setProblem(
    const Teuchos::ArrayView<const double>& source_centers,
    const Teuchos::ArrayView<const double>& target_centers,
    const double radius,
    const double alpha )
{
    DTK_REQUIRE( 0 == source_centers.size() % DIM );
    DTK_REQUIRE( 0 == target_centers.size() % DIM );

    // Build the local interpolation problems.
    buildLocalProblems( source_centers, target_centers, radius, alpha );

    DTK_ENSURE( Teuchos::nonnull(d_distributor) );
    DTK_ENSURE( d_local_problems.size() == target_centers.size() / DIM );
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
    const int target_lda,
    const int max_solve_iterations,
    const double solve_convergence_tolerance ) const
{
    DTK_REQUIRE( num_source_dims == num_target_dims );
    DTK_REQUIRE( source_data.size() == source_lda * num_source_dims );
    DTK_REQUIRE( target_data.size() == target_lda * num_target_dims );

    // Distribute the source data to the target decomposition.

    // Solve the local interpolation problems.
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the local interpolation problems.
 */
template<class Basis, class GO, int DIM>
void MovingLeastSquare<Basis,GO,DIM>::buildLocalProblems(
    const Teuchos::ArrayView<const double>& source_centers,
    const Teuchos::ArrayView<const double>& target_centers,
    const double radius,
    const double alpha )
{
    // Build the basis.
    Teuchos::RCP<Basis> basis = BP::create( radius );

    // Gather the source centers that are within a radius of the target
    // centers on this proc.
    Teuchos::Array<double> dist_sources;
    d_target_distributor = Teuchos::rcp( 
	new CenterDistributor<DIM>( 
	    d_comm, source_centers, target_centers, radius, dist_sources) );

    // Build the source/target pairings.
    SplineInterpolationPairing<DIM> target_pairings( 
	dist_sources, target_centers, radius );

    // Build the local problems.
    int num_targets = target_centers.size() / DIM;
    d_local_problems.resize( num_targets );
    Teuchos::ArrayView<const double> target_view;
    for ( int n = 0; n < num_targets; ++n )
    {
	target_view = target_centers(n*DIM,DIM);
	d_local_problems[n] = LocalMLSProblem<Basis,GO,DIM>(
	    target_view, target_pairings.childCenterIds(n),
	    dist_sources, *basis, alpha );
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MOVINGLEASTSQUARE_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_MovingLeastSquare_impl.hpp
//---------------------------------------------------------------------------//

