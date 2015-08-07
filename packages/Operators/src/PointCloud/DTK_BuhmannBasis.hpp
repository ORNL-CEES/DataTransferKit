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
 * \file   DTK_BuhmannBasis.hpp
 * \author Stuart R. Slattery
 * \brief  Buhmann compactly supported radial basis function.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BUHMANNBASIS_HPP
#define DTK_BUHMANNBASIS_HPP

#include "DTK_RadialBasisPolicy.hpp"

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class BuhmannBasis
 * \brief Buhmann compactly supported radial basis function.
 */
//---------------------------------------------------------------------------//
template<int ORDER>
class BuhmannBasis
{
  public:

    // Constructor.
    BuhmannBasis( const double radius )
	: d_radius( radius )
    { /* ... */ }

    // Compute the value of the basis at the given set of coordinates.
    double evaluateValue( const double x ) const;

    // Compute the gradient of the basis at the given set of coordinates.
    double evaluateGradient( const double x ) const;

  private:

    // Basis radius.
    double d_radius;
};

//---------------------------------------------------------------------------//
// RadialBasisPolicy implementation.
//---------------------------------------------------------------------------//
template<int ORDER>
class RadialBasisPolicy<BuhmannBasis<ORDER> >
{
  public:
    
    typedef BuhmannBasis<ORDER> spline_basis_type;
    
    static inline Teuchos::RCP<BuhmannBasis<ORDER> > 
    create( const double radius )
    { return Teuchos::rcp( new BuhmannBasis<ORDER>(radius) ); }

    static inline double evaluateValue( 
	const BuhmannBasis<ORDER>& basis, const double x )
    { return basis.evaluateValue( x ); }

    static inline double evaluateGradient( 
	const BuhmannBasis<ORDER>& basis, const double x )
    { return basis.evaluateGradient( x ); }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_BuhmannBasis_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_BUHMANNBASIS_HPP

//---------------------------------------------------------------------------//
// end DTK_BuhmannBasis.hpp
//---------------------------------------------------------------------------//

