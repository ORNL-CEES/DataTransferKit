/****************************************************************************
 * Copyright (c) 2012-2019 by the DataTransferKit authors                   *
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

// FIXME remove when able to use C++14
namespace std_ext
{
template <bool B, class T = void>
using enable_if_t = typename std::enable_if<B, T>::type;
}

//! Compute the maximum of two values.
template <typename T,
          typename = std_ext::enable_if_t<std::is_arithmetic<T>::value>>
KOKKOS_INLINE_FUNCTION T max( T a, T b )
{
    return ( a > b ) ? a : b;
}

//! Compute the minimum of two values.
template <typename T,
          typename = std_ext::enable_if_t<std::is_arithmetic<T>::value>>
KOKKOS_INLINE_FUNCTION T min( T a, T b )
{
    return ( a < b ) ? a : b;
}

/**
 * Branchless sign function. Return 1 if @param x is greater than zero, 0 if
 * @param x is zero, and -1 if @param x is less than zero.
 */
template <typename T,
          typename = std_ext::enable_if_t<std::is_arithmetic<T>::value>>
KOKKOS_INLINE_FUNCTION int sgn( T x )
{
    return ( x > 0 ) - ( x < 0 );
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

} // namespace KokkosHelpers
} // namespace DataTransferKit

#endif
