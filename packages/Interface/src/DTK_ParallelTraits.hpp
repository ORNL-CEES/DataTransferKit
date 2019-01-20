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
/*!
 * \file DTK_ParallelTraits.hpp
 * \brief Kokkos helpers.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_PARALLELTRAITS_HPP
#define DTK_PARALLELTRAITS_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_DefaultNode.hpp>

#include <type_traits>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Memory space aliases.

// Host memory.
#if defined( KOKKOS_ENABLE_SERIAL ) || defined( KOKKOS_ENABLE_OPENMP )
using HostSpace = Kokkos::HostSpace;
#endif

// Cuda Unified-Virtual-Memory
#if defined( KOKKOS_ENABLE_CUDA )
using CudaUVMSpace = Kokkos::CudaUVMSpace;
#endif

//---------------------------------------------------------------------------//
// Parallel execution space aliases.

// Serial
#if defined( KOKKOS_ENABLE_SERIAL )
using Serial = Kokkos::Serial;
#endif

// OpenMP
#if defined( KOKKOS_ENABLE_OPENMP )
using OpenMP = Kokkos::OpenMP;
#endif

// Cuda
#if defined( KOKKOS_ENABLE_CUDA )
using Cuda = Kokkos::Cuda;
#endif

//---------------------------------------------------------------------------//

} // namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_PARALLELTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_ParallelTraits.hpp
//---------------------------------------------------------------------------//
