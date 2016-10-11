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
 * \file   DTK_WuBasis_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Wu compactly supported radial basis function.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_WUBASIS_IMPL_HPP
#define DTK_WUBASIS_IMPL_HPP

#include <cmath>
#include <limits>

#include "DTK_EuclideanDistance.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Compute the value of the basis at the given set of
 * coordinates. Specialization of order 2.
 */
//---------------------------------------------------------------------------//
template <>
inline double WuBasis<2>::evaluateValue( const double radius,
                                         const double x ) const
{
    double xval = x / radius;
    double onemx = 1.0 - xval;
    double onemx2 = onemx * onemx;
    double xp2 = xval * xval;
    return ( xval < 1.0 )
               ? onemx * onemx2 * onemx2 *
                     ( 5.0 * xp2 * xp2 + 25.0 * xval * xp2 + 48.0 * xp2 +
                       40.0 * xval + 8.0 )
               : 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the gradient of the basis at the given set of
 * coordinates. Specialization of order 2.
 */
//---------------------------------------------------------------------------//
template <>
inline double WuBasis<2>::evaluateGradient( const double radius,
                                            const double x ) const
{
    double xval = x / radius;
    double xmone = xval - 1.0;
    double xmone2 = xmone * xmone;
    double xp2 = x * x;
    return ( xval < 1.0 )
               ? -9.0 * xmone2 * xmone2 * xval *
                     ( 5.0 * xp2 * xval + 20.0 * xp2 + 29.0 * xval + 16.0 )
               : 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the value of the basis at the given set of
 * coordinates. Specialization of order 4.
 */
//---------------------------------------------------------------------------//
template <>
inline double WuBasis<4>::evaluateValue( const double radius,
                                         const double x ) const
{
    double xval = x / radius;
    double onemx = 1.0 - xval;
    double onemx3 = onemx * onemx * onemx;
    double xp2 = xval * xval;
    double xp4 = xp2 * xp2;
    return ( xval < 1.0 )
               ? onemx3 * onemx3 *
                     ( 5.0 * xp4 * xval + 30.0 * xp4 + 72.0 * xp2 * xval +
                       82.0 * xp2 + 36.0 * xval + 6.0 )
               : 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the gradient of the basis at the given set of
 * coordinates. Specialization of order 4.
 */
//---------------------------------------------------------------------------//
template <>
inline double WuBasis<4>::evaluateGradient( const double radius,
                                            const double x ) const
{
    double xval = x / radius;
    double xmone = xval - 1.0;
    double xmone2 = xmone * xmone;
    double xp2 = xval * xval;
    return ( xval < 1.0 )
               ? 11.0 * xmone2 * xmone2 * xmone * xval *
                     ( 5.0 * xp2 * xp2 + 25.0 * xp2 * xval + 48.0 * xp2 +
                       40.0 * xval + 8.0 )
               : 0.0;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_WUBASIS_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_WuBasis_impl.hpp
//---------------------------------------------------------------------------//
