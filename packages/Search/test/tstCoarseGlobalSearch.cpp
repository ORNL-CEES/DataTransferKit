/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#include <Teuchos_UnitTestHarness.hpp>

#include <Teuchos_DefaultComm.hpp>

#include <DTK_CoarseGlobalSearch.hpp>
#include <DTK_DetailsPredicate.hpp>

namespace details = DataTransferKit::Details;

template <typename DeviceType>
class FillBoxes
{
  public:
    KOKKOS_INLINE_FUNCTION
    FillBoxes( Kokkos::View<DataTransferKit::Box *, DeviceType> boxes )
        : _boxes( boxes )
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( int const i ) const
    {
        if ( i == 0 )
            _boxes[0] = {0, 0, 0, 0, 0, 0};
        else
            _boxes[1] = {1, 1, 1, 1, 1, 1};
    }

  private:
    Kokkos::View<DataTransferKit::Box *, DeviceType> _boxes;
};

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CoarseGlobalSearch, constructor, NO )
{
    using namespace DataTransferKit;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();

    // Construct a set of local boxes.
    using DeviceType = typename CoarseGlobalSearch<NO>::DeviceType;
    using ExecutionSpace = typename DeviceType::execution_space;

    int const n = 2;
    Kokkos::View<Box *, DeviceType> boxes( "boxes", n );
    FillBoxes<DeviceType> fill_boxes_functor( boxes );
    Kokkos::parallel_for( "fill_boxes_functor",
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n ),
                          fill_boxes_functor );
    Kokkos::fence();

    CoarseGlobalSearch<NO> cgs( comm, boxes );
    (void)cgs;
}

// Include the test macros.
#include "DataTransferKitSearch_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CoarseGlobalSearch, constructor,     \
                                          NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
