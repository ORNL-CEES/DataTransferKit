/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef DTK_KOKKOS_HELPERS_HPP
#define DTK_KOKKOS_HELPERS_HPP

#include <Kokkos_Macros.hpp>

#include <cmath>   // isfinite
#include <cstdint> // uint32_t
#include <type_traits>

namespace DataTransferKit
{
namespace KokkosHelpers
{

// FIXME should be able to use directly Kokkos::ArithTrais<T>::infinity() when
// kokkos/kokkos-kernels#279 gets merged and eventually makes it into Trilinos.
template <typename T>
struct ArithTraits
{
    static KOKKOS_INLINE_FUNCTION T infinity();
};

template <>
struct ArithTraits<double>
{
    static KOKKOS_INLINE_FUNCTION double infinity() { return HUGE_VAL; }
};

//! Compute the maximum of two values.
template <typename T, typename = std::enable_if<std::is_arithmetic<T>::value>>
KOKKOS_INLINE_FUNCTION T max( T a, T b )
{
    return ( a > b ) ? a : b;
}

//! Compute the minimum of two values.
template <typename T, typename = std::enable_if<std::is_arithmetic<T>::value>>
KOKKOS_INLINE_FUNCTION T min( T a, T b )
{
    return ( a < b ) ? a : b;
}

/**
 * Branchless sign function. Return 1 if @param x is greater than zero, 0 if
 * @param x is zero, and -1 if @param x is less than zero.
 */
template <typename T, typename = std::enable_if<std::is_arithmetic<T>::value>>
KOKKOS_INLINE_FUNCTION int sgn( T x )
{
    return ( x > 0 ) - ( x < 0 );
}

/** Count the number of consecutive leading zero bits in 32 bit integer
 * @param x
 */
KOKKOS_INLINE_FUNCTION
int clz( uint32_t x )
{
#if defined( __CUDA_ARCH__ )
    // Note that the __clz() CUDA intrinsic function takes a signed integer
    // as input parameter.  This is fine but would need to be adjusted if
    // we were to change expandBits() and morton3D() to subdivide [0, 1]^3
    // into more 1024^3 bins.
    return __clz( x );
#elif defined( KOKKOS_COMPILER_GNU ) || ( KOKKOS_COMPILER_CLANG >= 500 )
    // According to https://en.wikipedia.org/wiki/Find_first_set
    // Clang 5.X supports the builtin function with the same syntax as GCC
    return ( x == 0 ) ? 32 : __builtin_clz( x );
#else
    if ( x == 0 )
        return 32;
    // The following is taken from:
    // http://stackoverflow.com/questions/23856596/counting-leading-zeros-in-a-32-bit-unsigned-integer-with-best-algorithm-in-c-pro
    const char debruijn32[32] = {0,  31, 9,  30, 3,  8,  13, 29, 2,  5, 7,
                                 21, 12, 24, 28, 19, 1,  10, 4,  14, 6, 22,
                                 25, 20, 11, 15, 23, 26, 16, 27, 17, 18};
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x++;
    return debruijn32[x * 0x076be629 >> 27];
#endif
}

/** Determine whether the given floating point argument @param x has finite
 * value.
 *
 * NOTE: Clang issues a warning if the std:: namespace is missing and nvcc
 * complains about calling a __host__ function from a __host__ __device__
 * function when it is present.
 */
template <typename FloatingPoint>
KOKKOS_INLINE_FUNCTION bool isFinite( FloatingPoint x )
{
#ifdef __CUDA_ARCH__
    return isfinite( x );
#else
    return std::isfinite( x );
#endif
}

/** Determine whether the given floating point number argument @param x is a
 * not-a-number (NaN) value.
 */
template <typename FloatingPoint>
KOKKOS_INLINE_FUNCTION bool isNan( FloatingPoint x )
{
#ifdef __CUDA_ARCH__
    return isnan( x );
#else
    return std::isnan( x );
#endif
}

} // namespace KokkosHelpers
} // namespace DataTransferKit

#endif
