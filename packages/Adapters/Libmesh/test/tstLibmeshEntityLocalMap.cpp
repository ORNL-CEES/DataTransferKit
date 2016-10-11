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
 * \file tstLibmeshEntityLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief LibmeshEntityLocalMap unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_LibmeshEntity.hpp>
#include <DTK_LibmeshEntityExtraData.hpp>
#include <DTK_LibmeshEntityLocalMap.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <libmesh/cell_hex8.h>
#include <libmesh/equation_systems.h>
#include <libmesh/libmesh.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/node.h>
#include <libmesh/parallel.h>
#include <libmesh/point.h>
#include <libmesh/system.h>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template <class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal>> getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp( new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// TEST EPSILON
//---------------------------------------------------------------------------//

const double epsilon = 1.0e-14;

//---------------------------------------------------------------------------//
// Hex-8 test.
TEUCHOS_UNIT_TEST( LibmeshEntity, hex_8_test )
{
    // Extract the raw mpi communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm<int>();
    Teuchos::RCP<const Teuchos::MpiComm<int>> mpi_comm =
        Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>( comm );
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> opaque_comm =
        mpi_comm->getRawMpiComm();
    MPI_Comm raw_comm = ( *opaque_comm )();

    // Create the mesh.
    int space_dim = 3;
    const std::string argv_string = "unit_test";
    const char *argv_char = argv_string.c_str();
    libMesh::LibMeshInit libmesh_init( 1, &argv_char, raw_comm );
    TEST_ASSERT( libMesh::initialized() );
    TEST_EQUALITY( (int)libmesh_init.comm().rank(), comm->getRank() );
    Teuchos::RCP<libMesh::Mesh> mesh =
        Teuchos::rcp( new libMesh::Mesh( libmesh_init.comm(), space_dim ) );

    // Create the nodes.
    int rank = comm->getRank();
    Teuchos::Array<libMesh::Node *> nodes( 8 );
    double node_coords[3];
    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = -2.0;
    nodes[0] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 0,
        rank );

    node_coords[0] = 2.0;
    node_coords[1] = 0.0;
    node_coords[2] = -2.0;
    nodes[1] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 1,
        rank );

    node_coords[0] = 2.0;
    node_coords[1] = 2.0;
    node_coords[2] = -2.0;
    nodes[2] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 2,
        rank );

    node_coords[0] = 0.0;
    node_coords[1] = 2.0;
    node_coords[2] = -2.0;
    nodes[3] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 3,
        rank );

    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;
    nodes[4] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 4,
        rank );

    node_coords[0] = 2.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;
    nodes[5] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 5,
        rank );

    node_coords[0] = 2.0;
    node_coords[1] = 2.0;
    node_coords[2] = 0.0;
    nodes[6] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 6,
        rank );

    node_coords[0] = 0.0;
    node_coords[1] = 2.0;
    node_coords[2] = 0.0;
    nodes[7] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 7,
        rank );

    // Make a hex-8.
    libMesh::Elem *hex_elem = mesh->add_elem( new libMesh::Hex8 );
    hex_elem->processor_id() = rank;
    hex_elem->set_id() = 2 * rank;
    for ( int i = 0; i < 8; ++i )
        hex_elem->set_node( i ) = nodes[i];

    // Check libmesh validity.
    mesh->libmesh_assert_valid_parallel_ids();

    // Make an adjacency data structure.
    DataTransferKit::LibmeshAdjacencies adjacencies( mesh );

    // Create a DTK entity for the hex.
    DataTransferKit::Entity dtk_entity =
        DataTransferKit::LibmeshEntity<libMesh::Elem>(
            Teuchos::ptr( hex_elem ), mesh.ptr(),
            Teuchos::ptrFromRef( adjacencies ) );

    // Make a libmesh system. We will put a first order linear basis on the
    // elements.
    libMesh::EquationSystems equation_systems( *mesh );
    libMesh::LinearImplicitSystem &system =
        equation_systems.add_system<libMesh::LinearImplicitSystem>( "Test" );
    system.add_variable( "test_var", libMesh::FIRST );

    // Create a local map from the libmesh mesh.
    Teuchos::RCP<DataTransferKit::EntityLocalMap> local_map =
        Teuchos::rcp( new DataTransferKit::LibmeshEntityLocalMap(
            mesh, Teuchos::rcpFromRef( system ) ) );

    // Test the measure.
    TEST_FLOATING_EQUALITY( local_map->measure( dtk_entity ), 8.0, epsilon );

    // Test the centroid.
    Teuchos::Array<double> centroid( space_dim, 0.0 );
    local_map->centroid( dtk_entity, centroid() );
    TEST_EQUALITY( centroid[0], 1.0 );
    TEST_EQUALITY( centroid[1], 1.0 );
    TEST_EQUALITY( centroid[2], -1.0 );

    // Make a good point and a bad point.
    Teuchos::Array<double> good_point( space_dim );
    good_point[0] = 0.5;
    good_point[1] = 1.5;
    good_point[2] = -1.0;
    Teuchos::Array<double> bad_point( space_dim );
    bad_point[0] = 0.75;
    bad_point[1] = -1.75;
    bad_point[2] = 0.35;

    // Test the reference frame safeguard.
    TEST_ASSERT(
        local_map->isSafeToMapToReferenceFrame( dtk_entity, good_point() ) );
    TEST_ASSERT(
        !local_map->isSafeToMapToReferenceFrame( dtk_entity, bad_point() ) );

    // Test the mapping to reference frame.
    Teuchos::Array<double> ref_good_point( space_dim );
    bool good_map = local_map->mapToReferenceFrame( dtk_entity, good_point(),
                                                    ref_good_point() );
    TEST_ASSERT( good_map );
    TEST_FLOATING_EQUALITY( ref_good_point[0], -0.5, epsilon );
    TEST_FLOATING_EQUALITY( ref_good_point[1], 0.5, epsilon );
    TEST_ASSERT( std::abs( ref_good_point[2] ) < epsilon );

    Teuchos::Array<double> ref_bad_point( space_dim );
    local_map->mapToReferenceFrame( dtk_entity, bad_point(), ref_bad_point() );

    // Test the point inclusion.
    TEST_ASSERT(
        local_map->checkPointInclusion( dtk_entity, ref_good_point() ) );
    TEST_ASSERT(
        !local_map->checkPointInclusion( dtk_entity, ref_bad_point() ) );

    // Test the map to physical frame.
    Teuchos::Array<double> phy_good_point( space_dim );
    local_map->mapToPhysicalFrame( dtk_entity, ref_good_point(),
                                   phy_good_point() );
    TEST_FLOATING_EQUALITY( good_point[0], phy_good_point[0], epsilon );
    TEST_FLOATING_EQUALITY( good_point[1], phy_good_point[1], epsilon );
    TEST_FLOATING_EQUALITY( good_point[2], phy_good_point[2], epsilon );

    Teuchos::Array<double> phy_bad_point( space_dim );
    local_map->mapToPhysicalFrame( dtk_entity, ref_bad_point(),
                                   phy_bad_point() );
    TEST_FLOATING_EQUALITY( bad_point[0], phy_bad_point[0], epsilon );
    TEST_FLOATING_EQUALITY( bad_point[1], phy_bad_point[1], epsilon );
    TEST_FLOATING_EQUALITY( bad_point[2], phy_bad_point[2], epsilon );

    // Test the coordinates of the points extracted through the centroid
    // function.
    DataTransferKit::Entity dtk_node;
    Teuchos::Array<double> point_coords( space_dim );
    int num_nodes = 8;
    for ( int n = 0; n < num_nodes; ++n )
    {
        dtk_node = DataTransferKit::LibmeshEntity<libMesh::Node>(
            Teuchos::ptr( nodes[n] ), mesh.ptr(),
            Teuchos::ptrFromRef( adjacencies ) );
        local_map->centroid( dtk_node, point_coords() );
        TEST_EQUALITY( ( *nodes[n] )( 0 ), point_coords[0] );
        TEST_EQUALITY( ( *nodes[n] )( 1 ), point_coords[1] );
        TEST_EQUALITY( ( *nodes[n] )( 2 ), point_coords[2] );
    }
}

//---------------------------------------------------------------------------//
// end tstLibmeshEntityLocalMap.cpp
//---------------------------------------------------------------------------//
