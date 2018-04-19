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
// Parallel node type aliases.

// Serial
#if defined( KOKKOS_HAVE_SERIAL )
using Serial = Kokkos::Serial;
#endif

// OpenMP
#if defined( KOKKOS_HAVE_OPENMP )
using OpenMP = Kokkos::OpenMP;
#endif

// Cuda
#if defined( KOKKOS_HAVE_CUDA )
using Cuda = Kokkos::Cuda;
#endif

//---------------------------------------------------------------------------//
// Check device-type to ensure Kokkos and Tpetra devices are compatible. This
// will not compile if they are incompatible.
template <class Device1, class Device2>
struct CompatibleDeviceTypes
{ /* ... */
};

template <class Device>
struct CompatibleDeviceTypes<Device, Device>
{
    using IsCompatible = std::true_type;
};

//---------------------------------------------------------------------------//
// Undefined Traits.
template <class Node>
class ParallelTraits
{ /* ... */
};

//---------------------------------------------------------------------------//
// Serial specialization.
#if defined( KOKKOS_HAVE_SERIAL )
template <>
class ParallelTraits<Serial>
{
  public:
    using ExecutionSpace = Kokkos::Serial;
    using DeviceType = typename ExecutionSpace::device_type;
    using MemorySpace = typename DeviceType::memory_space;
    using TpetraNode = ::Kokkos::Compat::KokkosSerialWrapperNode;

  private:
    using IsCompatible =
        typename CompatibleDeviceTypes<typename TpetraNode::device_type,
                                       DeviceType>::IsCompatible;
};
#endif

//---------------------------------------------------------------------------//
// OpenMP specialization.
#if defined( KOKKOS_HAVE_OPENMP )
template <>
class ParallelTraits<OpenMP>
{
  public:
    using ExecutionSpace = Kokkos::OpenMP;
    using DeviceType = typename ExecutionSpace::device_type;
    using MemorySpace = typename DeviceType::memory_space;
    using TpetraNode = ::Kokkos::Compat::KokkosOpenMPWrapperNode;

  private:
    using IsCompatible =
        typename CompatibleDeviceTypes<typename TpetraNode::device_type,
                                       DeviceType>::IsCompatible;
};
#endif

//---------------------------------------------------------------------------//
// Cuda specialization.
#if defined( KOKKOS_HAVE_CUDA )
template <>
class ParallelTraits<Cuda>
{
  public:
    using ExecutionSpace = Kokkos::Cuda;
    using DeviceType = typename ExecutionSpace::device_type;
    using MemorySpace = typename DeviceType::memory_space;
    using TpetraNode = ::Kokkos::Compat::KokkosCudaWrapperNode;

  private:
    using IsCompatible =
        typename CompatibleDeviceTypes<typename TpetraNode::device_type,
                                       DeviceType>::IsCompatible;
};
#endif

//---------------------------------------------------------------------------//

} // namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_PARALLELTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_ParallelTraits.hpp
//---------------------------------------------------------------------------//
