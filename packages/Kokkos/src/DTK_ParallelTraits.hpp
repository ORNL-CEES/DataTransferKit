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
    using TpetraNode = typename ::Kokkos::Compat::KokkosSerialWrapperNode;

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
    using TpetraNode = typename ::Kokkos::Compat::KokkosOpenMPWrapperNode;

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
    using TpetraNode = typename ::Kokkos::Compat::KokkosCudaWrapperNode;

  private:
    using IsCompatible =
        typename CompatibleDeviceTypes<typename TpetraNode::device_type,
                                       DeviceType>::IsCompatible;
};
#endif

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_PARALLELTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_ParallelTraits.hpp
//---------------------------------------------------------------------------//
