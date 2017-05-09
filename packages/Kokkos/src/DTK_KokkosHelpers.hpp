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

    /**
     * Branchless sign function. Return 1 if @param x is greater than zero, 0 if
     * @param x is zero, and -1 if @param x is less than zero.
     */
    KOKKOS_INLINE_FUNCTION
    static int sgn( int x ) { return ( x > 0 ) - ( x < 0 ); }
};

/**
 * This functor is similar to std::iota.
 */
template <typename NO, typename SC = int>
class Iota
{
  public:
    using DeviceType = typename NO::device_type;

    Iota( Kokkos::View<SC *, DeviceType> indices,
          SC offset = static_cast<SC>( 0 ) )
        : _indices( indices )
        , _offset( offset )
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( int const i ) const
    {
        _indices[i] = static_cast<SC>( i ) + _offset;
    }

  private:
    Kokkos::View<SC *, DeviceType> _indices;
    SC _offset;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_KOKKOSHELPERS_HPP

//---------------------------------------------------------------------------//
// end DTK_KokkosHelpers.hpp
//---------------------------------------------------------------------------//
