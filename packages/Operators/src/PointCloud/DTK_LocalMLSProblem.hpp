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
 * \file   DTK_LocalMLSProblem.hpp
 * \author Stuart R. Slattery
 * \brief  Local moving least square problem.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_LOCALMLSPROBLEM_HPP
#define DTK_LOCALMLSPROBLEM_HPP

#include "DTK_RadialBasisPolicy.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class LocalMLSProblem
 * \brief Local moving least square problem about a single target center using
 * quadratic polynomials.
 */
//---------------------------------------------------------------------------//
template<class Basis,int DIM>
class LocalMLSProblem
{
  public:

    //@{
    //! Typedefs.
    typedef RadialBasisPolicy<Basis> BP;
    //@}

    // Default constructor.
    LocalMLSProblem()
    { /* ... */ }

    // Constructor.
    LocalMLSProblem( const Teuchos::ArrayView<const double>& target_center,
		     const Teuchos::ArrayView<const unsigned>& source_lids,
		     const Teuchos::ArrayView<const double>& source_centers,
		     const Basis& basis,
		     const double radius );

    // Get a view of the local shape function.
    Teuchos::ArrayView<const double> shapeFunction() const
    { return d_shape_function(); }

  private:
    
    // Get a polynomial coefficient.
    double polynomialCoefficient( 
	const int coeff, const Teuchos::ArrayView<const double>& center ) const;

    // Check if a matrix is full rank.
    bool isFullRank( const Teuchos::SerialDenseMatrix<int,double>& matrix ) const;

  private:

    // Moving least square shape function.
    Teuchos::Array<double> d_shape_function;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_LOCALMLSPROBLEM_HPP

//---------------------------------------------------------------------------//
// end DTK_LocalMLSProblem.hpp
//---------------------------------------------------------------------------//

