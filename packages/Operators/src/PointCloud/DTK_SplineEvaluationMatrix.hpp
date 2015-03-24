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
 * \file   DTK_SplineEvaluationMatrix.hpp
 * \author Stuart R. Slattery
 * \brief  Spline transformation operator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SPLINEEVALUATIONMATRIX_HPP
#define DTK_SPLINEEVALUATIONMATRIX_HPP

#include "DTK_RadialBasisPolicy.hpp"
#include "DTK_SplineInterpolationPairing.hpp"
#include "DTK_PolynomialMatrix.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_RCP.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Operator.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class SplineEvaluationMatrix
 * \brief Sparse spline transformation operator (the A matrix). A = N + Q
 */
//---------------------------------------------------------------------------//
template<class Basis,int DIM>
class SplineEvaluationMatrix
{
  public:

    //@{
    //! Typedefs.
    typedef RadialBasisPolicy<Basis> BP;
    //@}

    // Constructor.
    SplineEvaluationMatrix(
	const Teuchos::RCP<const Tpetra::Map<int,DofId> >& domain_map,
	const Teuchos::RCP<const Tpetra::Map<int,DofId> >& range_map,
	const Teuchos::ArrayView<const double>& target_centers,
	const Teuchos::ArrayView<const DofId>& target_center_gids,
	const Teuchos::ArrayView<const double>& dist_source_centers,
	const Teuchos::ArrayView<const DofId>& dist_source_center_gids,
	const SplineInterpolationPairing<DIM>& target_pairings,
	const Basis& basis );

    // Get the basis component.
    Teuchos::RCP<Tpetra::Operator<double,int,DofId> > getN()
    { return d_N; }

    // Get the polynomial component.
    Teuchos::RCP<Tpetra::Operator<double,int,DofId> > getQ()
    { return d_Q; }

  private:

    // The N matrix.
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,DofId> > d_N;

    // The Q matrix.
    Teuchos::RCP<PolynomialMatrix> d_Q;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_SplineEvaluationMatrix_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_SPLINEEVALUATIONMATRIX_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineEvaluationMatrix.hpp
//---------------------------------------------------------------------------//

