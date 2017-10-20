/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#include "DTK_ConfigDefs.hpp"
#include "DTK_Types.h"

#include "PointCloudProblemGenerator/ExodusAsciiGenerator.hpp"
#include "PointCloudProblemGenerator/PointCloudProblemGenerator.hpp"

#include <Kokkos_View.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <memory>
#include <string>

//---------------------------------------------------------------------------//
// Partition the grids one-to-one
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ExodusAsciiGenerator, one_to_one, Node )
{
    // Type aliases.
    using DeviceType = typename Node::device_type;

    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Create a problem generator.
    std::string src_coord_file = "fine_sphere_node_coords.dat";
    std::string src_connectivity_file = "fine_sphere_node_connectivity.dat";
    std::string tgt_coord_file = "coarse_sphere_node_coords.dat";
    std::string tgt_connectivity_file = "coarse_sphere_node_connectivity.dat";
    auto generator = std::make_shared<
        DataTransferKit::ExodusAsciiGenerator<DeviceType, DeviceType>>(
        comm, src_coord_file, src_connectivity_file, tgt_coord_file,
        tgt_connectivity_file );

    // Generate a problem.
    Kokkos::View<DataTransferKit::Coordinate **, Kokkos::LayoutLeft, DeviceType>
        src_coords;
    Kokkos::View<DataTransferKit::Coordinate **, Kokkos::LayoutLeft, DeviceType>
        tgt_coords;
    generator->createOneToOneProblem( src_coords, tgt_coords );

    // Print sizes.
    std::cout << comm->getRank() << " " << src_coords.extent( 0 ) << " "
              << src_coords.extent( 1 ) << " " << tgt_coords.extent( 0 ) << " "
              << tgt_coords.extent( 1 ) << std::endl;
}

//---------------------------------------------------------------------------//
// Test partition the grids with ghosting.
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ExodusAsciiGenerator, ghosting, Node )
{
    // Type aliases.
    using DeviceType = typename Node::device_type;

    // Get the communicator.
    auto comm = Teuchos::DefaultComm<int>::getComm();

    // Create a problem generator.
    std::string src_coord_file = "fine_sphere_node_coords.dat";
    std::string src_connectivity_file = "fine_sphere_connectivity.dat";
    std::string tgt_coord_file = "coarse_sphere_node_coords.dat";
    std::string tgt_connectivity_file = "coarse_sphere_connectivity.dat";
    auto generator = std::make_shared<
        DataTransferKit::ExodusAsciiGenerator<DeviceType, DeviceType>>(
        comm, src_coord_file, src_connectivity_file, tgt_coord_file,
        tgt_connectivity_file );

    // Generate a problem.
    Kokkos::View<DataTransferKit::Coordinate **, Kokkos::LayoutLeft, DeviceType>
        src_coords;
    Kokkos::View<DataTransferKit::GlobalOrdinal *, Kokkos::LayoutLeft,
                 DeviceType>
        src_gids;
    Kokkos::View<DataTransferKit::Coordinate **, Kokkos::LayoutLeft, DeviceType>
        tgt_coords;
    Kokkos::View<DataTransferKit::GlobalOrdinal *, Kokkos::LayoutLeft,
                 DeviceType>
        tgt_gids;
    generator->createGeneralProblem( src_coords, src_gids, tgt_coords,
                                     tgt_gids );

    // Print sizes.
    std::cout << comm->getRank() << " " << src_coords.extent( 0 ) << " "
              << src_gids.extent( 0 ) << " " << tgt_coords.extent( 0 ) << " "
              << tgt_gids.extent( 0 ) << std::endl;
}
//---------------------------------------------------------------------------//

// Include the test macros.
#include "DataTransferKitDiscretization_ETIHelperMacros.h"

// Create the test group
#define UNIT_TEST_GROUP( NODE )                                                \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ExodusAsciiGenerator, one_to_one,    \
                                          NODE )                               \
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ExodusAsciiGenerator, ghosting, NODE )

// Demangle the types
DTK_ETI_MANGLING_TYPEDEFS()

// Instantiate the tests
DTK_INSTANTIATE_N( UNIT_TEST_GROUP )
