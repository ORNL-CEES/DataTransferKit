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
 * \file   DTK_MovingLeastSquares.hpp
 * \author Stuart R. Slattery
 * \brief Parallel moving least square interpolator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOVINGLEASTSQUARE_HPP
#define DTK_MOVINGLEASTSQUARE_HPP

#include "DTK_MeshFreeInterpolator.hpp"
#include "DTK_RadialBasisPolicy.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class MovingLeastSquare
 * \brief Parallel moving least square interpolator.
 *
 * The MovingLeastSquare is the top-level driver for parallel interpolation
 * problems.
 */
//---------------------------------------------------------------------------//
template<class Basis, class GO, int DIM>
class MovingLeastSquare : public MeshFreeInterpolator
{
  public:

    //@{
    //! Typedefs.
    typedef RadialBasisPolicy<Basis> BP;
    typedef Tpetra::Operator<double,int,GO> OP;
    typedef Tpetra::MultiVector<double,int,GO> MV;
    //@}

    // Constructor.
    MovingLeastSquare( const Teuchos::RCP<const Teuchos::Comm<int> >& comm );

    //! Destructor.
    ~MovingLeastSquare()
    { /* ... */ }

    // Set the interpolation problem.
    void setProblem( const Teuchos::ArrayView<const double>& source_centers,
		     const Teuchos::ArrayView<const double>& target_centers,
		     const double radius,
		     const double alpha );

    // Given a set of scalar values at the given source centers in the source
    // decomposition, interpolate them onto the target centers in the target
    // decomposition.
    void interpolate( const Teuchos::ArrayView<const double>& source_data,
		      const int num_source_dims,
		      const int source_lda,
		      const Teuchos::ArrayView<double>& target_data,
		      const int num_target_dims,
		      const int target_lda,
		      const int max_solve_iterations,
		      const double solve_convergence_tolerance ) const;

  private:

    // Build the local interpolation problems.
    void buildLocalProblems(
	const Teuchos::ArrayView<const double>& source_centers,
	const Teuchos::ArrayView<const double>& target_centers,
	const double radius,
	const double alpha );

  private:

    // Parallel communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > d_comm;

    // Source-to-target distributor.
    Teuchos::RCP<CenterDistributor<DIM> > d_distributor;

    // Local interpolation problems.
    Teuchos::Array<LocalMLSProblem<Basis,GO,DIM> > d_local_problems;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_MovingLeastSquare_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_MOVINGLEASTSQUARE_HPP

//---------------------------------------------------------------------------//
// end DTK_MovingLeastSquare.hpp
//---------------------------------------------------------------------------//

