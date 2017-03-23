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
 * \file   DTK_RadialBasisPolicy.hpp
 * \author Stuart R. Slattery
 * \brief  Policy class for spline interpolation basis functions.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_RADIALBASISPOLICY_HPP
#define DTK_RADIALBASISPOLICY_HPP

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Dummy struct. If a type does not create a specialization this will
 * not compile.
 */
template <typename UndefinedRadialBasis>
struct UndefinedRadialBasisPolicy
{
    static inline UndefinedRadialBasis notDefined()
    {
        return UndefinedRadialBasis::this_type_is_missing_a_specialization();
    }
};

//---------------------------------------------------------------------------//
/*!
 * \class RadialBasisPolicy
 * \brief Traits/policy class for compactly supported
 * radial basis functions spline interpolation.
 */
//---------------------------------------------------------------------------//
template <typename RadialBasis>
class RadialBasisPolicy
{
  public:
    //! Typedef for RadialBasis.
    typedef RadialBasis spline_basis_type;

    //! Creation method.
    static inline Teuchos::RCP<RadialBasis> create()
    {
        UndefinedRadialBasisPolicy<RadialBasis>::notDefined();
        return Teuchos::null;
    }

    //! Compute the value of the basis at the given value.
    static inline double evaluateValue( const RadialBasis &basis,
                                        const double radius, const double x )
    {
        UndefinedRadialBasisPolicy<RadialBasis>::notDefined();
        return 0.0;
    }

    //! Compute the gradient of the basis at the given value.
    static inline double evaluateGradient( const RadialBasis &basis,
                                           const double radius, const double x )
    {
        UndefinedRadialBasisPolicy<RadialBasis>::notDefined();
        return 0.0;
    }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_RADIALBASISPOLICY_HPP

//---------------------------------------------------------------------------//
// end DTK_RadialBasisPolicy.hpp
//---------------------------------------------------------------------------//
