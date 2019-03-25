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

#include <cmath> // isfinite
#include <type_traits>

namespace DataTransferKit
{
namespace KokkosHelpers
{

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
