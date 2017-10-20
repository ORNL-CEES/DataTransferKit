/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#include "DTK_ConfigDefs.hpp"

#include <Kokkos_View.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

//---------------------------------------------------------------------------//
// Test building and manipulating a tpetra vector with a kokkos view.
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraVectorView, vector_view, Node )
{
    // Template aliases.
    using SC = double;
    using LO = int;
    using GO = unsigned long int;
    using NO = Node;

    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Set the number of global elements.
    LO num_local = 100;
    GO num_global = num_local * comm->getSize();

    // Create a map.
    auto map =
        Tpetra::createUniformContigMapWithNode<LO, GO, NO>( num_global, comm );

    // Type aliases.
    using DualViewType =
        typename Tpetra::MultiVector<SC, LO, GO, NO>::dual_view_type;
    using DeviceType = typename DualViewType::device_type;
    using ExecutionSpace = typename DualViewType::execution_space;

    // Create a view and get the device view.
    int num_vec = 3;
    DualViewType dual_view( "test_dual_view", num_local, num_vec );
    auto dev_view = dual_view.template view<DeviceType>();

    // Put some data in the view.
    SC value = 3.3;
    auto fill_func = KOKKOS_LAMBDA( int i )
    {
        for ( int d = 0; d < num_vec; ++d )
            dev_view( i, d ) = value;
    };
    Kokkos::parallel_for( REGION_NAME( "fill_test_dual_vew" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, num_local ),
                          fill_func );

    // Create a vector from the view.
    auto vec = Teuchos::rcp(
        new Tpetra::MultiVector<SC, LO, GO, NO>( map, dual_view ) );

    // Check the vector norm.
    Kokkos::View<SC *, Kokkos::HostSpace> norms( "vector norms", num_vec );
    SC test_norm = num_global * value;
    vec->norm1( norms );
    for ( int d = 0; d < num_vec; ++d )
        TEST_FLOATING_EQUALITY( test_norm, norms( d ), 1.0e-14 );
}

//---------------------------------------------------------------------------//

// Include the test macros.
#include "DataTransferKitDiscretization_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraVectorView, vector_view, NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
