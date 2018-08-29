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
#ifndef DTK_PREDICATE_HPP
#define DTK_PREDICATE_HPP

#include <DTK_DetailsAlgorithms.hpp>
#include <DTK_DetailsNode.hpp>

namespace DataTransferKit
{
namespace Details
{
struct NearestPredicateTag
{
};
struct SpatialPredicateTag
{
};
} // namespace Details

template <typename Geometry>
struct Nearest
{
    using Tag = Details::NearestPredicateTag;

    KOKKOS_INLINE_FUNCTION
    Nearest() = default;

    KOKKOS_INLINE_FUNCTION
    Nearest( Geometry const &geometry, int k )
        : _geometry( geometry )
        , _k( k )
    {
    }

    Geometry _geometry;
    int _k = 0;
};

namespace Details
{

struct Sphere2
{
    KOKKOS_INLINE_FUNCTION
    Sphere2() = default;

    KOKKOS_INLINE_FUNCTION
    Sphere2( Point const &centroid, DistanceReturnType radius )
        : _centroid( centroid )
        , _radius( radius )
    {
    }

    KOKKOS_INLINE_FUNCTION
    Sphere2( Point const &centroid, double radius )
        : _centroid( centroid )
        , _radius( DistanceReturnType{radius * radius} )
    {
    }

    KOKKOS_INLINE_FUNCTION
    Point &centroid() { return _centroid; }

    KOKKOS_INLINE_FUNCTION
    Point const &centroid() const { return _centroid; }

    KOKKOS_INLINE_FUNCTION
    DistanceReturnType radius() const { return _radius; }

    Point _centroid;
    DistanceReturnType _radius = DistanceReturnType{0.};
};

KOKKOS_INLINE_FUNCTION
bool intersects( Sphere2 const &sphere, Box const &box )
{
    return distance( sphere.centroid(), box ) <= sphere.radius();
}

KOKKOS_INLINE_FUNCTION
Point return_centroid( Sphere2 const &sphere ) { return sphere.centroid(); }

KOKKOS_INLINE_FUNCTION
void expand( Box &box, Sphere2 const &sphere )
{
    for ( int d = 0; d < 3; ++d )
    {
        box.minCorner()[d] = KokkosHelpers::min(
            box.minCorner()[d], sphere.centroid()[d] - sphere.radius() );
        box.maxCorner()[d] = KokkosHelpers::max(
            box.maxCorner()[d], sphere.centroid()[d] + sphere.radius() );
    }
}

} // namespace Details

template <typename Geometry>
struct Intersects
{
    using Tag = Details::SpatialPredicateTag;

    KOKKOS_INLINE_FUNCTION Intersects() = default;

    KOKKOS_INLINE_FUNCTION Intersects( Geometry const &geometry )
        : _geometry( geometry )
    {
    }

    KOKKOS_INLINE_FUNCTION
    bool operator()( Node const *node ) const
    {
        return Details::intersects( _geometry, node->bounding_box );
    }

    Geometry _geometry;
};

using Within = Intersects<Details::Sphere2>;
using Overlap = Intersects<Box>;

template <typename Geometry>
KOKKOS_INLINE_FUNCTION Nearest<Geometry> nearest( Geometry const &geometry,
                                                  int k = 1 )
{
    return Nearest<Geometry>( geometry, k );
}

KOKKOS_INLINE_FUNCTION
Within within( Point const &p, double r )
{
    return Within( Details::Sphere2{p, Details::DistanceReturnType{r * r}} );
}

KOKKOS_INLINE_FUNCTION
Within within( Point const &p, Details::DistanceReturnType r )
{
    return Within( Details::Sphere2{p, r} );
}

KOKKOS_INLINE_FUNCTION
Overlap overlap( Box const &b ) { return Overlap( b ); }

} // namespace DataTransferKit

#endif
