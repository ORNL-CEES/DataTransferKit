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
 * \brief DTK_KokkosHelpers.hpp
 * \author Stuart R. Slattery
 * \brief Kokkos helpers.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_KOKKOSHELPERS_HPP
#define DTK_KOKKOSHELPERS_HPP

#include <Kokkos_Core.hpp>

#include <cfloat>
#include <climits>
#include <typeinfo>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class KokkosHelpers
  \brief Utility functions to help with Kokkos.
*/
//---------------------------------------------------------------------------//
class KokkosHelpers
{
  public:
    //! Compute the maximum of two values.
    template <class SC>
    KOKKOS_INLINE_FUNCTION static SC max( const SC left, const SC right )
    {
        return ( left > right ) ? left : right;
    }

    //! Compute the minimum of two values.
    template <class SC>
    KOKKOS_INLINE_FUNCTION static SC min( const SC left, const SC right )
    {
        return ( left < right ) ? left : right;
    }

    //! Return the maximum value representable by a given floating point or
    //! integer type. std::numeric_limits<SC>::max() is not available in CUDA
    //! 8 so we recreate that functionality here.
    template <class SC>
    KOKKOS_INLINE_FUNCTION static SC numericLimitsMax();

    //! Return the epsilon representable by a given floating point
    //! type. std::numeric_limits<SC>::epsilon() is not available in CUDA 8 so
    //! we recreate that functionality here.
    template <class SC>
    KOKKOS_INLINE_FUNCTION static SC numericLimitsEpsilon();
};

//---------------------------------------------------------------------------//
// Specializations of numericLimitsMax
//---------------------------------------------------------------------------//
// integer
template <>
KOKKOS_INLINE_FUNCTION int KokkosHelpers::numericLimitsMax<int>()
{
    return INT_MAX;
}

//---------------------------------------------------------------------------//
// unsigned integer
template <>
KOKKOS_INLINE_FUNCTION unsigned int
KokkosHelpers::numericLimitsMax<unsigned int>()
{
    return UINT_MAX;
}

//---------------------------------------------------------------------------//
// long integer
template <>
KOKKOS_INLINE_FUNCTION long int KokkosHelpers::numericLimitsMax<long int>()
{
    return LONG_MAX;
}

//---------------------------------------------------------------------------//
// unsigned long integer
template <>
KOKKOS_INLINE_FUNCTION unsigned long int
KokkosHelpers::numericLimitsMax<unsigned long int>()
{
    return ULONG_MAX;
}

//---------------------------------------------------------------------------//
// long long integer
template <>
KOKKOS_INLINE_FUNCTION long long int
KokkosHelpers::numericLimitsMax<long long int>()
{
    return LLONG_MAX;
}

//---------------------------------------------------------------------------//
// unsigned long long integer
template <>
KOKKOS_INLINE_FUNCTION unsigned long long int
KokkosHelpers::numericLimitsMax<unsigned long long int>()
{
    return ULLONG_MAX;
}

//---------------------------------------------------------------------------//
// float
template <>
KOKKOS_INLINE_FUNCTION float KokkosHelpers::numericLimitsMax<float>()
{
    return FLT_MAX;
}

//---------------------------------------------------------------------------//
// double
template <>
KOKKOS_INLINE_FUNCTION double KokkosHelpers::numericLimitsMax<double>()
{
    return DBL_MAX;
}

//---------------------------------------------------------------------------//
// long double
template <>
KOKKOS_INLINE_FUNCTION long double
KokkosHelpers::numericLimitsMax<long double>()
{
    return LDBL_MAX;
}

//---------------------------------------------------------------------------//
// Specializations of numericLimitsMax
//---------------------------------------------------------------------------//
// float
template <>
KOKKOS_INLINE_FUNCTION float KokkosHelpers::numericLimitsEpsilon<float>()
{
    return FLT_EPSILON;
}

//---------------------------------------------------------------------------//
// double
template <>
KOKKOS_INLINE_FUNCTION double KokkosHelpers::numericLimitsEpsilon<double>()
{
    return DBL_EPSILON;
}

//---------------------------------------------------------------------------//
// long double
template <>
KOKKOS_INLINE_FUNCTION long double
KokkosHelpers::numericLimitsEpsilon<long double>()
{
    return LDBL_EPSILON;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_KOKKOSHELPERS_HPP

//---------------------------------------------------------------------------//
// end DTK_KokkosHelpers.hpp
//---------------------------------------------------------------------------//
