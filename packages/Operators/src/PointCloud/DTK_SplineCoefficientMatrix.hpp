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
 * \file   DTK_SplineCoefficientMatrix.hpp
 * \author Stuart R. Slattery
 * \brief  Spline interpolation coefficient matrix.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SPLINECOEFFICIENTMATRIX_HPP
#define DTK_SPLINECOEFFICIENTMATRIX_HPP

#include "DTK_RadialBasisPolicy.hpp"
#include "DTK_SplineInterpolationPairing.hpp"
#include "DTK_PolynomialMatrix.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Operator.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class SplineCoefficientMatrix
 * \brief Sparse spline coefficient matrix.
 */
//---------------------------------------------------------------------------//
template<class Basis,int DIM>
class SplineCoefficientMatrix
{
  public:

    //@{
    //! Typedefs.
    typedef RadialBasisPolicy<Basis> BP;
    //@}

    // Constructor.
    SplineCoefficientMatrix(
	const Teuchos::RCP<const Tpetra::Map<int,SupportId> >& operator_map,
	const Teuchos::ArrayView<const double>& source_centers,
	const Teuchos::ArrayView<const SupportId>& source_center_gids,
	const Teuchos::ArrayView<const double>& dist_source_centers,
	const Teuchos::ArrayView<const SupportId>& dist_source_center_gids,
	const SplineInterpolationPairing<DIM>& source_pairings,
	const Basis& basis );

    // Get the basis component.
    Teuchos::RCP<Tpetra::Operator<double,int,SupportId> > getM()
    { return d_M; }

    // Get the polynomial component.
    Teuchos::RCP<Tpetra::Operator<double,int,SupportId> > getP()
    { return d_P; }

  private:

    // The M matrix.
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,SupportId> > d_M;

    // The P matrix.
    Teuchos::RCP<PolynomialMatrix> d_P;

};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_SplineCoefficientMatrix_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_SPLINECOEFFICIENTMATRIX_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineCoefficientMatrix.hpp
//---------------------------------------------------------------------------//

