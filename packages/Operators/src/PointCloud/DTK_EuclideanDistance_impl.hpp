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
 * \file   DTK_EuclideanDistance_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Wendland compactly supported radial basis function.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_EUCLIDEANDISTANCE_IMPL_HPP
#define DTK_EUCLIDEANDISTANCE_IMPL_HPP

#include <cmath>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Compute Euclidean distance between the given set of coordinates. 1-D
 * specialization.
 */
template <>
inline double EuclideanDistance<1>::distance( const double *x1,
                                              const double *x2 )
{
    return std::abs( x1[0] - x2[0] );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute Euclidean distance between the given set of coordinates. 2-D
 * specialization.
 */
template <>
inline double EuclideanDistance<2>::distance( const double *x1,
                                              const double *x2 )
{
    double xx = x1[0] - x2[0];
    double xy = x1[1] - x2[1];
    return std::sqrt( xx * xx + xy * xy );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute Euclidean distance between the given set of coordinates. 3-D
 * specialization.
 */
template <>
inline double EuclideanDistance<3>::distance( const double *x1,
                                              const double *x2 )
{
    double xx = x1[0] - x2[0];
    double xy = x1[1] - x2[1];
    double xz = x1[2] - x2[2];
    return std::sqrt( xx * xx + xy * xy + xz * xz );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_EUCLIDEANDISTANCE_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_EuclideanDistance_impl.hpp
//---------------------------------------------------------------------------//
