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
//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tstLibmeshEntityIntegrationRule.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  integration rule function test.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

#include "DTK_LibmeshAdjacencies.hpp"
#include "DTK_LibmeshEntity.hpp"
#include "DTK_LibmeshEntityIntegrationRule.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include <Teuchos_DefaultMpiComm.hpp>

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
// Hex-8 test.
TEUCHOS_UNIT_TEST( LibmeshEntityIntegrationRule, hex_8_test )
{
    // Extract the raw mpi communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
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
    node_coords[2] = 0.0;
    nodes[0] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 0,
        rank );

    node_coords[0] = 1.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;
    nodes[1] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 1,
        rank );

    node_coords[0] = 1.0;
    node_coords[1] = 1.0;
    node_coords[2] = 0.0;
    nodes[2] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 2,
        rank );

    node_coords[0] = 0.0;
    node_coords[1] = 1.0;
    node_coords[2] = 0.0;
    nodes[3] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 3,
        rank );

    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = 1.0;
    nodes[4] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 4,
        rank );

    node_coords[0] = 1.0;
    node_coords[1] = 0.0;
    node_coords[2] = 1.0;
    nodes[5] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 5,
        rank );

    node_coords[0] = 1.0;
    node_coords[1] = 1.0;
    node_coords[2] = 1.0;
    nodes[6] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 6,
        rank );

    node_coords[0] = 0.0;
    node_coords[1] = 1.0;
    node_coords[2] = 1.0;
    nodes[7] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ), 7,
        rank );

    // Make a hex-8.
    int num_nodes = 8;
    libMesh::Elem *hex_elem = mesh->add_elem( new libMesh::Hex8 );
    hex_elem->processor_id() = rank;
    hex_elem->set_id() = 2 * rank;
    for ( int i = 0; i < num_nodes; ++i )
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

    // Create an integration rule.
    Teuchos::RCP<DataTransferKit::EntityIntegrationRule> integration_rule =
        Teuchos::rcp( new DataTransferKit::LibmeshEntityIntegrationRule() );

    // Test the integration rule.
    Teuchos::Array<Teuchos::Array<double>> p_1;
    Teuchos::Array<double> w_1;
    integration_rule->getIntegrationRule( dtk_entity, 1, p_1, w_1 );
    TEST_EQUALITY( 1, w_1.size() );
    TEST_EQUALITY( 1, p_1.size() );
    TEST_EQUALITY( 3, p_1[0].size() );
    TEST_EQUALITY( 8.0, w_1[0] );
    TEST_EQUALITY( 0.0, p_1[0][0] );
    TEST_EQUALITY( 0.0, p_1[0][1] );
    TEST_EQUALITY( 0.0, p_1[0][2] );

    Teuchos::Array<Teuchos::Array<double>> p_2;
    Teuchos::Array<double> w_2;
    integration_rule->getIntegrationRule( dtk_entity, 2, p_2, w_2 );
    TEST_EQUALITY( 8, w_2.size() );
    TEST_EQUALITY( 8, p_2.size() );
    for ( int i = 0; i < 8; ++i )
    {
        TEST_EQUALITY( w_2[i], 1.0 );
        TEST_EQUALITY( p_2[i].size(), 3 );

        for ( int d = 0; d < 3; ++d )
        {
            TEST_FLOATING_EQUALITY( std::abs( p_2[i][d] ),
                                    1.0 / std::sqrt( 3.0 ), 1.0e-15 );
        }
    }
}

//---------------------------------------------------------------------------//
// end of tstLibmeshEntityIntegrationRule.cpp
//---------------------------------------------------------------------------//
