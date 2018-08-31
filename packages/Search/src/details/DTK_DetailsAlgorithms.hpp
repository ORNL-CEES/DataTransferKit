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
#ifndef DTK_DETAILS_ALGORITHMS_HPP
#define DTK_DETAILS_ALGORITHMS_HPP

#include <DTK_Box.hpp>
#include <DTK_KokkosHelpers.hpp> // isFinite, min, max
#include <DTK_Point.hpp>
#include <DTK_Sphere.hpp>

#include <Kokkos_Macros.hpp>

namespace DataTransferKit
{
namespace Details
{

KOKKOS_INLINE_FUNCTION
bool equals( Point const &l, Point const &r )
{
    for ( int d = 0; d < 3; ++d )
        if ( l[d] != r[d] )
            return false;
    return true;
}

KOKKOS_INLINE_FUNCTION
bool equals( Box const &l, Box const &r )
{
    return equals( l.minCorner(), r.minCorner() ) &&
           equals( l.maxCorner(), r.maxCorner() );
}

KOKKOS_INLINE_FUNCTION
bool equals( Sphere const &l, Sphere const &r )
{
    return equals( l.centroid(), r.centroid() ) && l.radius() == r.radius();
}

KOKKOS_INLINE_FUNCTION
bool isValid( Point const &p )
{
    using KokkosHelpers::isFinite;
    for ( int d = 0; d < 3; ++d )
        if ( !isFinite( p[d] ) )
            return false;
    return true;
}

KOKKOS_INLINE_FUNCTION
bool isValid( Box const &b )
{
    return isValid( b.minCorner() ) && isValid( b.maxCorner() );
}

KOKKOS_INLINE_FUNCTION
bool isValid( Sphere const &s )
{
    using KokkosHelpers::isFinite;
    return isValid( s.centroid() ) && isFinite( s.radius() ) &&
           ( s.radius() >= 0. );
}

// distance point-point
KOKKOS_INLINE_FUNCTION
double distance( Point const &a, Point const &b )
{
    double distance_squared = 0.0;
    for ( int d = 0; d < 3; ++d )
    {
        double tmp = b[d] - a[d];
        distance_squared += tmp * tmp;
    }
    return std::sqrt( distance_squared );
}

// distance point-box
KOKKOS_INLINE_FUNCTION
double distance( Point const &point, Box const &box )
{
    double s = 0;
    for ( int d = 0; d < 3; ++d )
    {
        double e = KokkosHelpers::max( box.minCorner()[d] - point[d], 0. ) +
                   KokkosHelpers::max( point[d] - box.maxCorner()[d], 0. );
        s += e * e;
    }
    return std::sqrt( s );
}

// distance point-sphere
KOKKOS_INLINE_FUNCTION
double distance( Point const &point, Sphere const &sphere )
{
    return KokkosHelpers::max(
        distance( point, sphere.centroid() ) - sphere.radius(), 0. );
}

// expand an axis-aligned bounding box to include a point
KOKKOS_INLINE_FUNCTION
void expand( Box &box, Point const &point )
{
    using KokkosHelpers::max;
    using KokkosHelpers::min;
    for ( int d = 0; d < 3; ++d )
    {
        box.minCorner()[d] = min( box.minCorner()[d], point[d] );
        box.maxCorner()[d] = max( box.maxCorner()[d], point[d] );
    }
}

// expand an axis-aligned bounding box to include another box
// NOTE: Box type is templated here to be able to use expand(box, box) in a
// Kokkos::parallel_reduce() in which case the arguments must be declared
// volatile.
template <typename BOX,
          typename = std::enable_if<std::is_same<
              typename std::remove_volatile<BOX>::type, Box>::value>>
KOKKOS_INLINE_FUNCTION void expand( BOX &box, BOX const &other )
{
    for ( int d = 0; d < 3; ++d )
    {
        box.minCorner()[d] =
            KokkosHelpers::min( box.minCorner()[d], other.minCorner()[d] );
        box.maxCorner()[d] =
            KokkosHelpers::max( box.maxCorner()[d], other.maxCorner()[d] );
    }
}

// expand an axis-aligned bounding box to include a sphere
KOKKOS_INLINE_FUNCTION
void expand( Box &box, Sphere const &sphere )
{
    for ( int d = 0; d < 3; ++d )
    {
        box.minCorner()[d] = KokkosHelpers::min(
            box.minCorner()[d], sphere.centroid()[d] - sphere.radius() );
        box.maxCorner()[d] = KokkosHelpers::max(
            box.maxCorner()[d], sphere.centroid()[d] + sphere.radius() );
    }
}

// check if two axis-aligned bounding boxes intersect
KOKKOS_INLINE_FUNCTION
bool intersects( Box const &box, Box const &other )
{
    double s = 0;
    for ( int d = 0; d < 3; ++d )
        s +=
            KokkosHelpers::max( box.minCorner()[d] - other.maxCorner()[d],
                                0. ) +
            KokkosHelpers::max( other.minCorner()[d] - box.maxCorner()[d], 0. );
    return ( s == 0 );
}

// check if a sphere intersects with an  axis-aligned bounding box
KOKKOS_INLINE_FUNCTION
bool intersects( Sphere const &sphere, Box const &box )
{
    Point const &c = sphere.centroid();
    double r = sphere.radius();
    double s = 0.;
    for ( int d = 0; d < 3; ++d )
    {
        double e = KokkosHelpers::max( box.minCorner()[d] - c[d], 0. ) +
                   KokkosHelpers::max( c[d] - box.maxCorner()[d], 0. );
        s += e * e;
    }
#if defined( __CUDA_ARCH__ )
    // WTF for CUDA sqrt is faster this way????
    // return ( s <= r * r );               // slow
    // return ( std::sqrt( s ) <= r );      // faster
    return ( s / r <= r ); // the fastest
#else
    // FIXME: this breaks DistributedSearchTree tests due to inconsistency with
    // distance() which does take the sqrt. But this is much faster
    return ( s <= r * r );
#endif
}

// calculate the centroid of a box
KOKKOS_INLINE_FUNCTION
void centroid( Box const &box, Point &c )
{
    for ( int d = 0; d < 3; ++d )
        c[d] = 0.5 * ( box.minCorner()[d] + box.maxCorner()[d] );
}

KOKKOS_INLINE_FUNCTION
Point return_centroid( Point const &point ) { return point; }

KOKKOS_INLINE_FUNCTION
Point return_centroid( Box const &box )
{
    Point c;
    for ( int d = 0; d < 3; ++d )
        c[d] = 0.5 * ( box.minCorner()[d] + box.maxCorner()[d] );
    return c;
}

KOKKOS_INLINE_FUNCTION
Point return_centroid( Sphere const &sphere ) { return sphere.centroid(); }

} // namespace Details
} // namespace DataTransferKit

#endif
