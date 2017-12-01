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

    // Create a view and get the device view.
    int num_vec = 3;
    DualViewType dual_view( "test_dual_view", num_local, num_vec );
    auto dev_view = dual_view.template view<DeviceType>();

    // Put some data in the view.
    SC value = 3.3;
    Kokkos::deep_copy( dev_view, value );

    // Create a vector from the view.
    auto vec = Teuchos::rcp(
        new Tpetra::MultiVector<SC, LO, GO, NO>( map, dual_view ) );

    // Check the vector norm.
    Kokkos::View<SC *, DeviceType> norms( "vector norms", num_vec );
    SC test_norm = num_global * value;
    vec->norm1( norms );
    Kokkos::View<SC *, Kokkos::HostSpace> host_norms( "host vector norms",
                                                      num_vec );
    Kokkos::deep_copy( host_norms, norms );
    for ( int d = 0; d < num_vec; ++d )
        TEST_FLOATING_EQUALITY( test_norm, host_norms( d ), 1.0e-14 );
}

//---------------------------------------------------------------------------//
// Test converting between a regular view and a dual view.
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TpetraVectorView, dual_view, Node )
{
    // Tpetra template aliases.
    using SC = double;
    using LO = int;
    using GO = unsigned long int;
    using NO = Node;

    // Type aliases.
    using DeviceType = typename Node::device_type;
    using DualViewType =
        typename Tpetra::MultiVector<SC, LO, GO, NO>::dual_view_type;
    using HostType = typename DualViewType::t_host::device_type;

    // Create a dual view compatible with a Tpetra vector from the regular
    // view on the device. Start by creating a host view compatible with the
    // device view. Note the create_mirror_view is being called here instead
    // of create_mirror. By calling create_mirror_view if the host and device
    // memory spaces are the same but the execution spaces are different it
    // will simply return the original view. In the case that they are
    // different it will allocate a new view of the device (e.g. the CUDA
    // case). This is the only way this will work with the DualView as the
    // DualView sync semantics only do the deep copy if the memory spaces are
    // different.
    auto createDualViewFromDeviceView =
        []( const Kokkos::View<SC **, Kokkos::LayoutLeft, DeviceType>
                &d_view ) {
            auto h_view = Kokkos::create_mirror_view( d_view );
            return DualViewType( d_view, h_view );
        };

    // Get a device view from a dual view.
    auto getDeviceViewFromDualView = []( const DualViewType &dual_view ) {
        return dual_view.template view<DeviceType>();
    };

    // ---------------------------------------------------------------------------
    // PART 1: Create a device view, and create a dual view from that device
    // view.

    // Create a device view.
    LO num_local = 10;
    int num_vec = 3;
    Kokkos::View<SC **, Kokkos::LayoutLeft, DeviceType> device_view(
        "device_view", num_local, num_vec );

    // Then create the dual view.
    DualViewType dual_view = createDualViewFromDeviceView( device_view );

    // -----------------------------------------------------------------------
    // PART 2: Fill the device view with data, sync the host and device, and
    // make sure that the host data changed to be the same as the device data.

    // Put some data in the regular view on the device.
    SC value = 3.3;
    Kokkos::deep_copy( device_view, value );

    // Mark the device data in the dual view as modified.
    dual_view.template modify<DeviceType>();

    // Sync the dual view - this will copy the data from the device to the
    // host. This only happens if the host and device memory spaces are
    // actually different. If they are the same then the view created by
    // create_mirror_view points to the same data as the device view and
    // therefore a copy is not necessary. If you created the host view with
    // create_mirror instead this test would fail for non-GPU cases because
    // that view would be full of zeros from initialization and would not get
    // a call to deep_copy because the memory spaces are the same.
    dual_view.template sync<HostType>();

    // Get the host view out of the device.
    Kokkos::View<SC **, Kokkos::LayoutLeft, DeviceType> host_view =
        dual_view.template view<HostType>();

    // Check that the host data was assigned the scalar value.
    for ( int i = 0; i < num_local; ++i )
        for ( int d = 0; d < num_vec; ++d )
            TEST_EQUALITY( value, host_view( i, d ) );

    // ----------------------------------------------------------------------
    // PART 3: Change the data in the dual view on the host and sync the host
    // and device.

    // Now change the host data.
    SC value_2 = 2.12;
    Kokkos::deep_copy( host_view, value_2 );

    // Mark the host data in the dual view as modified.
    dual_view.template modify<HostType>();

    // Copy the data from the host to the device.
    dual_view.template sync<DeviceType>();

    // ---------------------------------------------------------------------
    // PART 4: Create a device view from the dual view, copy the data in the
    // device view to a host view, and check that the device data was synced
    // with the new data.

    // Get a device view from the dual view.
    Kokkos::View<SC **, Kokkos::LayoutLeft, DeviceType> device_view_2 =
        getDeviceViewFromDualView( dual_view );

    // Copy the device view to a fresh host view.
    Kokkos::View<SC **, Kokkos::LayoutLeft, HostType> host_view_2(
        "host_view_t", num_local, num_vec );
    Kokkos::deep_copy( host_view_2, device_view_2 );

    // Check that the host data was assigned the new_scalar value.
    for ( int i = 0; i < num_local; ++i )
        for ( int d = 0; d < num_vec; ++d )
            TEST_EQUALITY( value_2, host_view_2( i, d ) );
}

//---------------------------------------------------------------------------//

// Include the test macros.
#include "DataTransferKitUtils_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraVectorView, vector_view,       \
                                          NODE )                               \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( TpetraVectorView, dual_view, NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
