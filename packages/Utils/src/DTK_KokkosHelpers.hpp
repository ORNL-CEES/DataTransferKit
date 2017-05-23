/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/
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

    /** Count the number of consecutive leading zero bits in 32 bit integer
     * @param x.
     */
    KOKKOS_INLINE_FUNCTION
    static int clz( uint32_t x )
    {
        if ( x == 0 )
            return 32;
        // The following is taken from:
        // http://stackoverflow.com/questions/23856596/counting-leading-zeros-in-a-32-bit-unsigned-integer-with-best-algorithm-in-c-pro
        static const char debruijn32[32] = {
            0, 31, 9, 30, 3, 8,  13, 29, 2,  5,  7,  21, 12, 24, 28, 19,
            1, 10, 4, 14, 6, 22, 25, 20, 11, 15, 23, 26, 16, 27, 17, 18};
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        x++;
        return debruijn32[x * 0x076be629 >> 27];
    }
};

/**
 * This functor is similar to std::iota.
 */
template <typename DeviceType, typename SC = int>
class Iota
{
  public:
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
