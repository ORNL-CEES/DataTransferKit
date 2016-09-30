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
 * \file   DTK_BuhmannBasis_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Buhmann compactly supported radial basis function.
 */
//---------------------------------------------------------------------------//

#ifndef DTKBUHMANNBASIS_IMPL_HPP
#define DTKBUHMANNBASIS_IMPL_HPP

#include <cmath>
#include <limits>

#include "DTK_EuclideanDistance.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Compute the value of the basis at the given set of
 * coordinates. Specialization of order 3.
 */
//---------------------------------------------------------------------------//
template<>
inline double
BuhmannBasis<3>::evaluateValue( const double radius, const double x ) const
{
    double xval = x / radius;
    double xp2 = xval*xval;
    double xp4 = xp2*xp2;
    double xp6 = xp4*xp2;
    double xp72 = std::sqrt(xp6*xval);
    double q1 = 84.0/5.0;
    double q2 = 1024.0/5.0;
    return ( x > radius ) ? 0.0 :
        xp4*xp4 - q1*xp6 + q2*xval*xp72 - 378.0*xp4 + q2*xp72 - q1*xp2 + 1.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the gradient of the basis at the given set of
 * coordinates. Specialization of order 3.
 */
//---------------------------------------------------------------------------//
template<>
inline double
BuhmannBasis<3>::evaluateGradient( const double radius, const double x ) const
{
    double xval = x / radius;
    double xp2 = xval*xval;
    double xp3 = xp2*xval;
    double xp5 = xp3*xp2;
    double xp52 = std::sqrt(xp5);
    return ( x > radius ) ? 0.0 :
        (8.0/5.0) * (576.0*xval*xp52 + 448.0*xp52 + 5.0*xp5*xp2 -
                     63.0*xp5 - 945.0*xp3 - 21.0*xval);
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTKBUHMANNBASIS_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_BuhmannBasis_impl.hpp
//---------------------------------------------------------------------------//

