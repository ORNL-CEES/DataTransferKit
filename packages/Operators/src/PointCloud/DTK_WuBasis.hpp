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
 * \file   DTK_WuBasis.hpp
 * \author Stuart R. Slattery
 * \brief  Wu compactly supported radial basis function.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_WUBASIS_HPP
#define DTK_WUBASIS_HPP

#include "DTK_RadialBasisPolicy.hpp"

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class WuBasis
 * \brief Wu compactly supported radial basis function.
 */
//---------------------------------------------------------------------------//
template<int ORDER>
class WuBasis
{
  public:

    // Compute the value of the basis at the given set of coordinates.
    double evaluateValue( const double radius, const double x ) const;

    // Compute the gradient of the basis at the given set of coordinates.
    double evaluateGradient( const double radius, const double x ) const;
};

//---------------------------------------------------------------------------//
// RadialBasisPolicy implementation.
//---------------------------------------------------------------------------//
template<int ORDER>
class RadialBasisPolicy<WuBasis<ORDER> >
{
  public:

    typedef WuBasis<ORDER> spline_basis_type;

    static inline Teuchos::RCP<WuBasis<ORDER> > create()
    { return Teuchos::rcp( new WuBasis<ORDER>() ); }

    static inline double evaluateValue(
        const WuBasis<ORDER>& basis, const double radius, const double x )
    { return basis.evaluateValue( radius, x ); }

    static inline double evaluateGradient(
        const WuBasis<ORDER>& basis, const double radius, const double x )
    { return basis.evaluateGradient( radius, x ); }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_WuBasis_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_WUBASIS_HPP

//---------------------------------------------------------------------------//
// end DTK_WuBasis.hpp
//---------------------------------------------------------------------------//

