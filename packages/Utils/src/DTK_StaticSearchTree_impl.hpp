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
 * \file DTK_StaticSearchTree_impl.hpp
 * \author Stuart R. Slattery
 * \brief Spatial searching for point clouds.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_STATICSEARCHTREE_IMPL_HPP
#define DTK_STATICSEARCHTREE_IMPL_HPP

#include <limits>

#include "DTK_DBC.hpp"

#include <Teuchos_as.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Compute the distance between a given point and a point in the cloud.
 */
template <>
inline double
PointCloud<1>::kdtree_distance( const double *p1, const std::size_t idx_p2,
                                std::size_t DTK_REMEMBER( size ) ) const
{
    DTK_REQUIRE( 1 == size );
    DTK_REQUIRE( idx_p2 < Teuchos::as<std::size_t>( d_points.size() ) );
    const double d0 = p1[0] - d_points[idx_p2];
    return d0 * d0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the distance between a given point and a point in the cloud.
 */
template <>
inline double
PointCloud<2>::kdtree_distance( const double *p1, const std::size_t idx_p2,
                                std::size_t DTK_REMEMBER( size ) ) const
{
    DTK_REQUIRE( 2 == size );
    DTK_REQUIRE( 2 * idx_p2 + 1 < Teuchos::as<std::size_t>( d_points.size() ) );
    const double d0 = p1[0] - d_points[2 * idx_p2];
    const double d1 = p1[1] - d_points[2 * idx_p2 + 1];
    return d0 * d0 + d1 * d1;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the distance between a given point and a point in the cloud.
 */
template <>
inline double
PointCloud<3>::kdtree_distance( const double *p1, const std::size_t idx_p2,
                                std::size_t DTK_REMEMBER( size ) ) const
{
    DTK_REQUIRE( 3 == size );
    DTK_REQUIRE( 3 * idx_p2 + 2 < Teuchos::as<std::size_t>( d_points.size() ) );
    const double d0 = p1[0] - d_points[3 * idx_p2];
    const double d1 = p1[1] - d_points[3 * idx_p2 + 1];
    const double d2 = p1[2] - d_points[3 * idx_p2 + 2];
    return d0 * d0 + d1 * d1 + d2 * d2;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the point coordinate at the given dimension.
 */
template <>
inline double PointCloud<1>::kdtree_get_pt( const std::size_t idx,
                                            int DTK_REMEMBER( dim ) ) const
{
    DTK_REQUIRE( dim < 1 );
    DTK_REQUIRE( idx < Teuchos::as<std::size_t>( d_points.size() ) );
    return d_points[idx];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the point coordinate at the given dimension.
 */
template <>
inline double PointCloud<2>::kdtree_get_pt( const std::size_t idx,
                                            int dim ) const
{
    DTK_REQUIRE( dim < 2 );
    DTK_REQUIRE( 2 * idx + dim < Teuchos::as<std::size_t>( d_points.size() ) );
    return d_points[2 * idx + dim];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the point coordinate at the given dimension.
 */
template <>
inline double PointCloud<3>::kdtree_get_pt( const std::size_t idx,
                                            int dim ) const
{
    DTK_REQUIRE( dim < 3 );
    DTK_REQUIRE( 3 * idx + dim < Teuchos::as<std::size_t>( d_points.size() ) );
    return d_points[3 * idx + dim];
}

//---------------------------------------------------------------------------//
// NanoflannTree Implementation.
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template <int DIM>
NanoflannTree<DIM>::NanoflannTree(
    const Teuchos::ArrayView<const double> &points,
    const unsigned max_leaf_size )
{
    DTK_CHECK( 0 == points.size() % DIM );

    d_cloud = PointCloud<DIM>( points );
    d_tree = Teuchos::rcp( new TreeType(
        DIM, d_cloud,
        nanoflann::KDTreeSingleIndexAdaptorParams( max_leaf_size ) ) );
    d_tree->buildIndex();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform an n-nearest neighbor search.
 */
template <int DIM>
Teuchos::Array<unsigned>
NanoflannTree<DIM>::nnSearch( const Teuchos::ArrayView<const double> &point,
                              const unsigned num_neighbors ) const
{
    DTK_REQUIRE( DIM == point.size() );
    Teuchos::Array<unsigned> neighbors( num_neighbors );
    Teuchos::Array<double> neighbor_dists( num_neighbors );
    d_tree->knnSearch( point.getRawPtr(), num_neighbors, neighbors.getRawPtr(),
                       neighbor_dists.getRawPtr() );
    return neighbors;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Perform a nearest neighbor search within a specified radius.
 */
template <int DIM>
Teuchos::Array<unsigned>
NanoflannTree<DIM>::radiusSearch( const Teuchos::ArrayView<const double> &point,
                                  const double radius ) const
{
    DTK_REQUIRE( DIM == point.size() );
    Teuchos::Array<std::pair<unsigned, double>> neighbor_pairs;
    nanoflann::SearchParams params;
    double l2_radius = radius * radius;
    d_tree->radiusSearch( point.getRawPtr(), l2_radius, neighbor_pairs,
                          params );

    Teuchos::Array<std::pair<unsigned, double>>::const_iterator pair_it;
    Teuchos::Array<unsigned> neighbors( neighbor_pairs.size() );
    Teuchos::Array<unsigned>::iterator id_it;
    for ( id_it = neighbors.begin(), pair_it = neighbor_pairs.begin();
          id_it != neighbors.end(); ++id_it, ++pair_it )
    {
        *id_it = pair_it->first;
    }

    return neighbors;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_STATICSEARCHTREE_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_StaticSearchTree_impl.hpp
//---------------------------------------------------------------------------//
