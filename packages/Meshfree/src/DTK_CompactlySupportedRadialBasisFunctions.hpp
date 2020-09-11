/****************************************************************************
 * Copyright (c) 2012-2020 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef DTK_COMPACTLY_SUPPORTED_RADIAL_BASIS_FUNCTIONS_HPP
#define DTK_COMPACTLY_SUPPORTED_RADIAL_BASIS_FUNCTIONS_HPP

#include <Kokkos_Macros.hpp>

#include <cmath> // log, sqrt

namespace DataTransferKit
{

template <typename RBF>
class RadialBasisFunction
{
  public:
    KOKKOS_FUNCTION
    RadialBasisFunction( double radius )
        : _radius( radius )
    {
        // FIXME check precondition radius is greater than zero
    }
    KOKKOS_INLINE_FUNCTION double operator()( double x ) const
    {
        return _rbf( x / _radius );
    }

  private:
    double _radius;
    RBF _rbf;
};

template <int k>
struct Wendland;

template <>
struct Wendland<0>
{
    KOKKOS_INLINE_FUNCTION double operator()( double x ) const
    {
        return ( 1.0 - x ) * ( 1.0 - x );
    }
};

template <>
struct Wendland<2>
{
    KOKKOS_INLINE_FUNCTION double operator()( double x ) const
    {
        return ( 1.0 - x ) * ( 1.0 - x ) * ( 1.0 - x ) * ( 1.0 - x ) *
               ( 4.0 * x + 1.0 );
    }
};

template <>
struct Wendland<4>
{
    KOKKOS_INLINE_FUNCTION double operator()( double x ) const
    {
        return ( 1.0 - x ) * ( 1.0 - x ) * ( 1.0 - x ) * ( 1.0 - x ) *
               ( 1.0 - x ) * ( 1.0 - x ) * ( 35.0 * x * x + 18.0 * x + 3.0 );
    }
};

template <>
struct Wendland<6>
{
    KOKKOS_INLINE_FUNCTION double operator()( double x ) const
    {
        return ( 1.0 - x ) * ( 1.0 - x ) * ( 1.0 - x ) * ( 1.0 - x ) *
               ( 1.0 - x ) * ( 1.0 - x ) * ( 1.0 - x ) * ( 1.0 - x ) *
               ( 32.0 * x * x * x + 25.0 * x * x + 8.0 * x + 1.0 );
    }
};

template <int k>
struct Wu;

template <>
struct Wu<2>
{
    KOKKOS_INLINE_FUNCTION double operator()( double x ) const
    {
        return ( 1.0 - x ) * ( 1.0 - x ) * ( 1.0 - x ) * ( 1.0 - x ) *
               ( 3.0 * x * x * x + 12.0 * x * x + 16.0 * x + 4.0 );
    }
};

template <>
struct Wu<4>
{
    KOKKOS_INLINE_FUNCTION double operator()( double x ) const
    {
        return ( 1.0 - x ) * ( 1.0 - x ) * ( 1.0 - x ) * ( 1.0 - x ) *
               ( 1.0 - x ) * ( 1.0 - x ) *
               ( 5.0 * x * x * x * x * x + 30.0 * x * x * x * x +
                 72.0 * x * x * x + 82.0 * x * x + 36.0 * x + 6.0 );
    }
};

template <int k>
struct Buhmann;

template <>
struct Buhmann<2>
{
    KOKKOS_INLINE_FUNCTION double operator()( double x ) const
    {
        return 2.0 * x * x * x * x * log( x ) - 7.0 / 2.0 * x * x * x * x +
               16 / 3.0 * x * x * x - 2 * x * x + 1.0 / 6.0;
    }
};

template <>
struct Buhmann<3>
{
    KOKKOS_INLINE_FUNCTION double operator()( double x ) const
    {
        return x * x * x * x * x * x * x * x -
               84.0 / 5.0 * x * x * x * x * x * x +
               1024.0 / 5.0 * x * x * x * x * sqrt( x ) -
               378.0 * x * x * x * x + 1024.0 / 5.0 * x * x * x * sqrt( x ) -
               84.0 / 5.0 * x * x + 1.0;
    }
};

template <>
struct Buhmann<4>
{
    KOKKOS_INLINE_FUNCTION double operator()( double x ) const
    {
        return 99.0 / 35.0 * x * x * x * x * x * x * x * x -
               132.0 * x * x * x * x * x * x +
               9216.0 / 35.0 * x * x * x * x * x * sqrt( x ) -
               11264.0 / 35.0 * x * x * x * x * sqrt( x ) +
               198.0 * x * x * x * x - 396.0 / 5.0 * x * x + 1.0;
    }
};

} // namespace DataTransferKit

#endif
