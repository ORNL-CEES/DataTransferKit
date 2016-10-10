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
 * \file   DTK_WendlandBasis_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Wendland compactly supported radial basis function.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_WENDLANDBASIS_IMPL_HPP
#define DTK_WENDLANDBASIS_IMPL_HPP

#include <cmath>
#include <limits>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Compute the value of the basis at the given set of
 * coordinates. Specialization of order 0.
 */
//---------------------------------------------------------------------------//
template<>
inline double
WendlandBasis<0>::evaluateValue( const double radius, const double x ) const
{
    double xval = x / radius;
    double onemx = 1.0 - xval;
    return ( xval < 1.0 ) ? onemx*onemx : 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the gradient of the basis at the given set of
 * coordinates. Specialization of order 0.
 */
//---------------------------------------------------------------------------//
template<>
inline double
WendlandBasis<0>::evaluateGradient( const double radius, const double x ) const
{
    double xval = x / radius;
    return ( xval < 1.0 ) ? 2.0*xval - 2.0 : 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the value of the basis at the given set of
 * coordinates. Specialization of order 2.
 */
//---------------------------------------------------------------------------//
template<>
inline double
WendlandBasis<2>::evaluateValue( const double radius, const double x ) const
{
    double xval = x / radius;
    double onemx = 1.0 - xval;
    double onemx2 = onemx*onemx;
    return ( xval < 1.0 ) ? onemx2*onemx2*(4.0*xval + 1.0) : 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the gradient of the basis at the given set of
 * coordinates. Specialization of order 2.
 */
//---------------------------------------------------------------------------//
template<>
inline double
WendlandBasis<2>::evaluateGradient( const double radius, const double x ) const
{
    double xval = x / radius;
    double xmone = xval - 1.0;
    return ( xval < 1.0 ) ? 20.0*xval*xmone*xmone*xmone : 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the value of the basis at the given set of
 * coordinates. Specialization of order 4.
 */
//---------------------------------------------------------------------------//
template<>
inline double
WendlandBasis<4>::evaluateValue( const double radius, const double x ) const
{
    double xval = x / radius;
    double onemx = 1.0 - xval;
    double onemx2 = onemx*onemx;
    return ( xval < 1.0 ) ?
        onemx2*onemx2*onemx2*(35.0*xval*xval + 18.0*xval + 3.0) : 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the gradient of the basis at the given set of
 * coordinates. Specialization of order 4.
 */
//---------------------------------------------------------------------------//
template<>
inline double
WendlandBasis<4>::evaluateGradient( const double radius, const double x ) const
{
    double xval = x / radius;
    double xmone = xval - 1.0;
    double xmone2 = xmone*xmone;
    return ( xval < 1.0 ) ? 56.0*xval*xmone2*xmone2*xmone*(5.0*xval+1) : 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the value of the basis at the given set of
 * coordinates. Specialization of order 6.
 */
//---------------------------------------------------------------------------//
template<>
inline double
WendlandBasis<6>::evaluateValue( const double radius, const double x ) const
{
    double xval = x / radius;
    double onemx = 1.0 - xval;
    double onemx2 = onemx*onemx;
    double onemx4 = onemx2*onemx2;
    double xs = xval*xval;
    return ( xval < 1.0 ) ?
        onemx4*onemx4*(32.0*xs*xval+25.0*xs+8.0*xval+1.0) : 0.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the gradient of the basis at the given set of
 * coordinates. Specialization of order 6.
 */
//---------------------------------------------------------------------------//
template<>
inline double
WendlandBasis<6>::evaluateGradient( const double radius, const double x ) const
{
    double xval = x / radius;
    double xmone = xval - 1.0;
    double xmone2 = xmone*xmone;
    double xmone4 = xmone2*xmone2;
    return ( xval < 1.0 ) ?
        22.0*xval*xmone4*xmone2*xmone*(16.0*xval*xval+7.0*xval+1.0) : 0.0;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_WENDLANDBASIS_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_WendlandBasis_impl.hpp
//---------------------------------------------------------------------------//

