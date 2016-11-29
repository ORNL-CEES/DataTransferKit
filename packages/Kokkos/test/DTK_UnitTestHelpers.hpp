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
 * \file   DTK_UnitTestHelpers.hpp
 * \author Stuart Slattery
 * \brief  Unit test helpers.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_UNITTESTHELPERS_HPP
#define DTK_UNITTESTHELPERS_HPP

#include <Kokkos_Core.hpp>

#include <Teuchos_UnitTestHarness.hpp>

//---------------------------------------------------------------------------//
// Instantiate the unit tests with different execution spaces. The input
// TEST_GROUP_FUNC is a unary function which takes an execution space
// argument. Developers should wrap their test template instantiations in such
// a unary function to automatically activate various node types.

// serial
#ifdef KOKKOS_HAVE_SERIAL
using serial_device = Kokkos::Serial;
#define DTK_SERIAL_TEST_INSTANT( TEST_GROUP_FUNC ) \
    TEST_GROUP_FUNC( serial_device )
#else
#define DTK_SERIAL_TEST_INSTANT( TEST_GROUP_FUNC )
#endif

// pthread
#ifdef KOKKOS_HAVE_PTHREAD
using pthread_device = Kokkos::Threads;
#define DTK_PTHREAD_TEST_INSTANT( TEST_GROUP_FUNC ) \
    TEST_GROUP_FUNC( pthread_device )
#else
#define DTK_PTHREAD_TEST_INSTANT( TEST_GROUP_FUNC )
#endif

// open mp
#ifdef KOKKOS_HAVE_OPENMP
using openmp_device = Kokkos::OpenMP;
#define DTK_OPENMP_TEST_INSTANT( TEST_GROUP_FUNC ) \
    TEST_GROUP_FUNC( openmp_device )
#else
#define DTK_OPENMP_TEST_INSTANT( TEST_GROUP_FUNC )
#endif

// cuda
#ifdef KOKKOS_HAVE_CUDA
using cuda_device = Kokkos::Cuda;
#define DTK_CUDA_TEST_INSTANT( TEST_GROUP_FUNC ) \
    TEST_GROUP_FUNC( cuda_device )
#else
#define DTK_CUDA_TEST_INSTANT( TEST_GROUP_FUNC )
#endif

// call all instantiations
#define DTK_TEST_DEVICE_INSTANT( TEST_GROUP_FUNC ) \
    DTK_SERIAL_TEST_INSTANT( TEST_GROUP_FUNC )          \
    DTK_PTHREAD_TEST_INSTANT( TEST_GROUP_FUNC )         \
    DTK_OPENMP_TEST_INSTANT( TEST_GROUP_FUNC )          \
    DTK_CUDA_TEST_INSTANT( TEST_GROUP_FUNC )

//---------------------------------------------------------------------------//
    
#endif // end DTK_UNITTESTHELPERS_HPP

//---------------------------------------------------------------------------//
// end DTK_UnitTestHelpers.hpp
//---------------------------------------------------------------------------//
