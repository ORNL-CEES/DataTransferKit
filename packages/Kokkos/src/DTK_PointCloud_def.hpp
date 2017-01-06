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
 * \brief DTK_PointCloud_def.hpp
 * \author Stuart R. Slattery
 * \brief Point cloud.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_POINTCLOUD_DEF_HPP
#define DTK_POINTCLOUD_DEF_HPP

#include "DTK_ConfigDefs.hpp"
#include "DTK_DBC.hpp"
#include "DTK_KokkosHelpers.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Ptr.hpp>

#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_Core.hpp>

#include <type_traits>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Local bounding box Kokkos functor.
//---------------------------------------------------------------------------//
//! Functor to compute the local bounding box of the point cloud.
template <class SC, class LO, class GO, class NO>
class LocalBoundingBoxFunctor
{
  public:
    //! Coordinate type alias.
    using coordinate_view =
        typename PointCloud<SC, LO, GO, NO>::coordinate_view;
    typedef SC value_type[];
    using size_type = typename coordinate_view::traits::size_type;

    //! Constructor.
    LocalBoundingBoxFunctor( const coordinate_view coordinates )
        : _coords( coordinates )
        , value_count( 6 )
    { /* ... */
    }

    //! Reduction operator.
    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type p, value_type update ) const
    {
        for ( size_t d = 0; d < _coords.extent( 1 ); ++d )
        {
            update[d] = KokkosHelpers::min( _coords( p, d ), update[d] );
            update[d + 3] =
                KokkosHelpers::max( _coords( p, d ), update[d + 3] );
        }
    }

    //! Join operator.
    KOKKOS_INLINE_FUNCTION
    void join( volatile value_type dst, const volatile value_type src ) const
    {
        for ( size_t d = 0; d < _coords.extent( 1 ); ++d )
        {
            dst[d] = KokkosHelpers::min( src[d], dst[d] );
            dst[d + 3] = KokkosHelpers::max( src[d + 3], dst[d + 3] );
        }
    }

    //! Init operator
    KOKKOS_INLINE_FUNCTION
    void init( value_type dst ) const
    {
        for ( size_t d = 0; d < 3; ++d )
        {
            dst[d] = Kokkos::Details::ArithTraits<SC>::max();
            dst[d + 3] = -Kokkos::Details::ArithTraits<SC>::max();
        }
    }

  public:
    // Point cloud coordinates.
    coordinate_view _coords;

    // Box size - specifically named this way for Kokkos
    size_type value_count;
};

//---------------------------------------------------------------------------//
// Class Implementation.
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 * \param comm Parallel communicator.
 * \param global_ids Point global ids. Dimensions: (Point)
 * \param coordinates Point coordinates. Dimensions: (Point,SpaceDim)
*/
template <class SC, class LO, class GO, class NO>
PointCloud<SC, LO, GO, NO>::PointCloud(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
    const global_id_view global_ids, const coordinate_view coordinates )
    : _comm( comm )
    , _global_ids( global_ids )
    , _coordinates( coordinates )
{
    DTK_REQUIRE( _global_ids.extent( 0 ) == _coordinates.extent( 0 ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the spatial dimension of the point cloud.
 * \return The spatial dimension of the point cloud.
 */
template <class SC, class LO, class GO, class NO>
size_t PointCloud<SC, LO, GO, NO>::spaceDim() const
{
    return _coordinates.extent( 1 );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the local number of points in the point cloud.
 * \return The local number of points in the point cloud.
 */
template <class SC, class LO, class GO, class NO>
size_t PointCloud<SC, LO, GO, NO>::numLocalPoints() const
{
    return _global_ids.extent( 0 );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global number of points in the point cloud.
 * \return The global number of points in the point cloud.
 */
template <class SC, class LO, class GO, class NO>
global_size_t PointCloud<SC, LO, GO, NO>::numGlobalPoints() const
{
    size_t local_num = numLocalPoints();
    size_t global_num = 0;
    Teuchos::reduceAll( *_comm, Teuchos::REDUCE_SUM, local_num,
                        Teuchos::ptrFromRef( global_num ) );
    return global_num;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the local bounding box.
 * \return The local bounding box (min x, min y, min z, max x, max y, max z).
 */
template <class SC, class LO, class GO, class NO>
Kokkos::Array<SC, 6> PointCloud<SC, LO, GO, NO>::localBoundingBox() const
{
    SC local_box[6];

    LocalBoundingBoxFunctor<SC, LO, GO, NO> functor( _coordinates );

    Kokkos::parallel_reduce(
        "Point cloud local bounding box",
        Kokkos::RangePolicy<execution_space>( 0, numLocalPoints() ), functor,
        local_box );
    Kokkos::fence();

    Kokkos::Array<SC, 6> box_array;
    for ( int i = 0; i < 6; ++i )
        box_array[i] = local_box[i];
    return box_array;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global bounding box.
 * \return The global bounding box (min x, min y, min z, max x, max y, max z).
 */
template <class SC, class LO, class GO, class NO>
Kokkos::Array<SC, 6> PointCloud<SC, LO, GO, NO>::globalBoundingBox() const
{
    Kokkos::Array<SC, 6> local_box = localBoundingBox();
    Kokkos::Array<SC, 6> global_box;

    for ( int i = 3; i < 6; ++i )
        local_box[i] *= -1;

    Teuchos::reduceAll( *_comm, Teuchos::REDUCE_MIN, 6, &local_box[0],
                        &global_box[0] );

    for ( int i = 3; i < 6; ++i )
        global_box[i] *= -1;

    return global_box;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the point global ids.
 * \return The point global ids. Dimensions: (Point)
 */
template <class SC, class LO, class GO, class NO>
auto PointCloud<SC, LO, GO, NO>::globalIds() const -> const global_id_view
{
    return _global_ids;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the point coordinates.
 * \return The point coordinates. Dimensions: (Point,SpaceDim)
 */
template <class SC, class LO, class GO, class NO>
auto PointCloud<SC, LO, GO, NO>::coordinates() const -> const coordinate_view
{
    return _coordinates;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Explicit Instantiation Macro.
//---------------------------------------------------------------------------//
#define DTK_POINTCLOUD_INSTANT( SCALAR, LO, GO, NODE )                         \
    template class PointCloud<SCALAR, LO, GO, NODE>;

//---------------------------------------------------------------------------//

#endif // end DTK_POINTCLOUD_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_PointCloud_def.hpp
//---------------------------------------------------------------------------//
