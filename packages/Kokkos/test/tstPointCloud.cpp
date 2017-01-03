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
 * \file   tstPointCloud.cpp
 * \author Stuart Slattery
 * \brief  Point Cloud unit tests.
 */
//---------------------------------------------------------------------------//

#include <DTK_PointCloud.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include <limits>

//---------------------------------------------------------------------------//
// FUNCTOR HELPERS
//---------------------------------------------------------------------------//
//! Create global ids.
template <class SC, class LO, class GO, class NO>
class CreateIdsFunctor
{
  public:
    using PointCloud = DataTransferKit::PointCloud<SC, LO, GO, NO>;
    using DeviceType = typename NO::device_type;
    using GlobalIdView = Kokkos::View<GO *, DeviceType>;

    CreateIdsFunctor( GlobalIdView global_ids, int comm_rank )
        : _global_ids( global_ids )
        , _comm_rank( comm_rank )
    { /* ... */
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_t p ) const
    {
        _global_ids( p ) = _comm_rank * _global_ids.extent( 0 ) + p;
    }

  private:
    GlobalIdView _global_ids;
    int _comm_rank;
};

//---------------------------------------------------------------------------//
// TEST TEMPLATE DECLARATIONS
//---------------------------------------------------------------------------//
// Point cloud test.
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( PointCloud, basic, SC, LO, GO, NO )
{
    // Test types.
    using PointCloud = DataTransferKit::PointCloud<SC, LO, GO, NO>;
    using DeviceType = typename NO::device_type;
    using GlobalIdView = Kokkos::View<GO *, DeviceType>;
    using CoordinateView = Kokkos::View<SC **, DeviceType>;
    using ExecutionSpace = typename PointCloud::execution_space;

    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Create global ids.
    size_t local_num_points = 10000;
    GlobalIdView global_ids( "global_ids", local_num_points );
    CreateIdsFunctor<SC, LO, GO, NO> id_functor( global_ids, comm_rank );
    Kokkos::parallel_for(
        "fill_gids", Kokkos::RangePolicy<ExecutionSpace>( 0, local_num_points ),
        id_functor );
    Kokkos::fence();

    // Create point coordinates.
    size_t space_dim = 3;
    CoordinateView coordinates( "coordinates", local_num_points, space_dim );
    Kokkos::Random_XorShift64_Pool<ExecutionSpace> rng_pool( 34998321 );
    for ( size_t d = 0; d < 3; ++d )
    {
        auto sv = Kokkos::subview( coordinates, Kokkos::ALL(), d );
        Kokkos::fill_random( sv, rng_pool, local_num_points );
        Kokkos::fence();
    }

    // Create a point cloud.
    PointCloud cloud( comm, global_ids, coordinates );

    // Test the point cloud data.
    TEST_EQUALITY( cloud.spaceDim(), space_dim );
    TEST_EQUALITY( cloud.numLocalPoints(), local_num_points );
    TEST_EQUALITY( cloud.numGlobalPoints(), comm_size * local_num_points );

    // Check the global ids.
    auto cloud_ids = cloud.globalIds();
    auto host_ids = Kokkos::create_mirror_view( cloud_ids );
    Kokkos::deep_copy( host_ids, cloud_ids );
    for ( size_t p = 0; p < local_num_points; ++p )
    {
        typename GlobalIdView::value_type test_id =
            comm_rank * local_num_points + p;
        TEST_EQUALITY( host_ids( p ), test_id );
    }

    // Compute the local bounding box manually on the host.
    auto cloud_coords = cloud.coordinates();
    auto host_coords = Kokkos::create_mirror_view( cloud_coords );
    Kokkos::deep_copy( host_coords, cloud_coords );
    Kokkos::Array<SC, 6> host_local_box;
    for ( size_t d = 0; d < space_dim; ++d )
    {
        host_local_box[d] = std::numeric_limits<SC>::max();
        host_local_box[d + 3] = -std::numeric_limits<SC>::max();
    }
    for ( size_t p = 0; p < local_num_points; ++p )
    {
        for ( size_t d = 0; d < space_dim; ++d )
        {
            host_local_box[d] =
                std::min( host_local_box[d], host_coords( p, d ) );
            host_local_box[d + 3] =
                std::max( host_local_box[d + 3], host_coords( p, d ) );
        }
    }

    // Compute the global bounding box manually on the host.
    Kokkos::Array<SC, 6> host_global_box;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_MIN, 3, &host_local_box[0],
                        &host_global_box[0] );

    Teuchos::reduceAll( *comm, Teuchos::REDUCE_MAX, 3, &host_local_box[3],
                        &host_global_box[3] );

    // Check the bounding boxes.
    auto local_box = cloud.localBoundingBox();
    TEST_COMPARE_ARRAYS( local_box, host_local_box );

    auto global_box = cloud.globalBoundingBox();
    TEST_COMPARE_ARRAYS( global_box, host_global_box );
}

//---------------------------------------------------------------------------//
// TEST TEMPLATE INSTANTIATIONS
//---------------------------------------------------------------------------//

// Include the test macros.
#include "DataTransferKitKokkos_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE )                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( PointCloud, basic, SCALAR, LO, GO,   \
                                          NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

//---------------------------------------------------------------------------//
// end tstPointCloud.cpp
//---------------------------------------------------------------------------//
