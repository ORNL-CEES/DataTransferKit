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
 * \file   DTK_NewtonSolver_impl.hpp
 * \author Stuart Slattery
 * \brief  A stateless class for Newton's method.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_NEWTONSOLVER_IMPL_HPP
#define DTK_NEWTONSOLVER_IMPL_HPP

#include <cmath>
#include <limits>

#include "DTK_DBC.hpp"

#include <Teuchos_ScalarTraits.hpp>

#include <Intrepid_RealSpaceTools.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Get the center of the reference cell of the given topology.
 */
template<typename NonlinearProblem>
void NewtonSolver<NonlinearProblem>::solve( MDArray& u, 
					    NonlinearProblem& problem,
					    const double tolerance,
					    const int max_iters )
{
    DTK_REQUIRE( 2 == u.rank() );

    // Allocate nonlinear residual, Jacobian, Newton update, and work arrays.
    int d0 = u.dimension(0);
    int d1 = u.dimension(1);
    MDArray F( d0, d1 );
    MDArray J( d0, d1, d1 );
    MDArray J_inv( d0, d1, d1 );
    MDArray update( d0, d1 );
    MDArray u_old = u;
    MDArray conv_check( d0 );

    // Compute the initial state.
    NPT::updateState( problem, u );

    // Computen the initial nonlinear residual and scale by -1 to get -F(u).
    NPT::evaluateResidual( problem, u, F );
    Intrepid::RealSpaceTools<Scalar>::scale( 
	F, -Teuchos::ScalarTraits<Scalar>::one() );

    // Compute the initial Jacobian.
    NPT::evaluateJacobian( problem, u, J );

    // Check for degeneracy of the Jacobian. If it is degenerate then the
    // problem is ill conditioned and return very large numbers in the state
    // vector that correspond to no solution.
    MDArray det( 1 );
    Intrepid::RealSpaceTools<Scalar>::det( det, J );
    if ( std::abs(det(0)) < tolerance )
    {
	for ( int m = 0; m < d0; ++m )
	{
	    for ( int n = 0; n < d1; ++n )
	    {
		u(m,n) = std::numeric_limits<Scalar>::max();
	    }
	}
	return;
    }

    // Nonlinear solve.
    for ( int k = 0; k < max_iters; ++k )
    {
	// Solve the linear model, delta_u = J^-1 * -F(u).
	Intrepid::RealSpaceTools<Scalar>::inverse( J_inv, J );
	Intrepid::RealSpaceTools<Scalar>::matvec( update, J_inv, F );

	// Update the solution, u += delta_u.
	Intrepid::RealSpaceTools<Scalar>::add( u, update );

	// Check for convergence.
	Intrepid::RealSpaceTools<Scalar>::subtract( u_old, u );
	Intrepid::RealSpaceTools<Scalar>::vectorNorm( 
	    conv_check, u_old, Intrepid::NORM_TWO );
	if ( tolerance > conv_check(0) )
	{
	    break;
	}

	// Reset for the next iteration.
	u_old = u;

	// Update any state-dependent data from the last iteration using the
	// new solution vector.
	NPT::updateState( problem, u );

	// Compute the nonlinear residual and scale by -1 to get -F(u).
	NPT::evaluateResidual( problem, u, F );
	Intrepid::RealSpaceTools<Scalar>::scale( 
	    F, -Teuchos::ScalarTraits<Scalar>::one() );

	// Compute the Jacobian.
	NPT::evaluateJacobian( problem, u, J );
    }

    // Check for convergence.
    DTK_ENSURE( tolerance > conv_check(0) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_NEWTONSOLVER_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_NewtonSolver_impl.hpp
//---------------------------------------------------------------------------//

