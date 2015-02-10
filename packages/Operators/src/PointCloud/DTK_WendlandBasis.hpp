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
 * \file   DTK_WendlandBasis.hpp
 * \author Stuart R. Slattery
 * \brief  Wendland compactly supported radial basis function.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_WENDLANDBASIS_HPP
#define DTK_WENDLANDBASIS_HPP

#include "DTK_RadialBasisPolicy.hpp"

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class WendlandBasis
 * \brief Wendland compactly supported radial basis function.
 */
//---------------------------------------------------------------------------//
template<int ORDER>
class WendlandBasis
{
  public:

    // Constructor.
    WendlandBasis( const double radius )
	: d_radius( radius )
    { /* ... */ }

    // Destructor.
    ~WendlandBasis()
    { /* ... */ }

    // Compute the value of the basis at the given value.
    double evaluateValue( const double x ) const;

    // Compute the gradient of the basis at the given value.
    double evaluateGradient( const double x ) const;

  private:

    // Basis radius.
    double d_radius;
};

//---------------------------------------------------------------------------//
// RadialBasisPolicy implementation.
//---------------------------------------------------------------------------//
template<int ORDER>
class RadialBasisPolicy<WendlandBasis<ORDER> >
{
  public:
    
    typedef WendlandBasis<ORDER> spline_basis_type;
    
    static inline Teuchos::RCP<WendlandBasis<ORDER> > 
    create( const double radius )
    { return Teuchos::rcp( new WendlandBasis<ORDER>(radius) ); }

    static inline double evaluateValue( 
	const WendlandBasis<ORDER>& basis, const double x )
    { return basis.evaluateValue( x ); }

    static inline double evaluateGradient( 
	const WendlandBasis<ORDER>& basis, const double x )
    { return basis.evaluateGradient( x ); }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_WendlandBasis_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_WENDLANDBASIS_HPP

//---------------------------------------------------------------------------//
// end DTK_WendlandBasis.hpp
//---------------------------------------------------------------------------//

