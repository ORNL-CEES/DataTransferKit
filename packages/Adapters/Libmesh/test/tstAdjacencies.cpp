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
 * \file tstLibmeshAdjacencies.cpp
 * \author Stuart R. Slattery
 * \brief LibmeshAdjacencies unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_LibmeshAdjacencies.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <libmesh/cell_hex8.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/node.h>
#include <libmesh/parallel.h>
#include <libmesh/point.h>

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
// Hex-8 test.
TEUCHOS_UNIT_TEST( LibmeshAdjacencies, hex_8_test )
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
    const std::string argv_string = "--keep-cout";
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
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ),
        8 * rank + 0, rank );

    node_coords[0] = 1.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;
    nodes[1] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ),
        8 * rank + 1, rank );

    node_coords[0] = 1.0;
    node_coords[1] = 1.0;
    node_coords[2] = 0.0;
    nodes[2] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ),
        8 * rank + 2, rank );

    node_coords[0] = 0.0;
    node_coords[1] = 1.0;
    node_coords[2] = 0.0;
    nodes[3] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ),
        8 * rank + 3, rank );

    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = 1.0;
    nodes[4] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ),
        8 * rank + 4, rank );

    node_coords[0] = 1.0;
    node_coords[1] = 0.0;
    node_coords[2] = 1.0;
    nodes[5] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ),
        8 * rank + 5, rank );

    node_coords[0] = 1.0;
    node_coords[1] = 1.0;
    node_coords[2] = 1.0;
    nodes[6] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ),
        8 * rank + 6, rank );

    node_coords[0] = 0.0;
    node_coords[1] = 1.0;
    node_coords[2] = 1.0;
    nodes[7] = mesh->add_point(
        libMesh::Point( node_coords[0], node_coords[1], node_coords[2] ),
        8 * rank + 7, rank );

    // Make a hex-8.
    libMesh::Elem *hex_elem_1 = mesh->add_elem( new libMesh::Hex8 );
    hex_elem_1->processor_id() = rank;
    hex_elem_1->set_id() = 2 * rank;
    for ( int i = 0; i < 8; ++i )
        hex_elem_1->set_node( i ) = nodes[i];

    // Make another hex-8.
    libMesh::Elem *hex_elem_2 = mesh->add_elem( new libMesh::Hex8 );
    hex_elem_2->processor_id() = rank;
    hex_elem_2->set_id() = 2 * rank + 1;
    for ( int i = 0; i < 8; ++i )
        hex_elem_2->set_node( i ) = nodes[i];

    // Check libmesh validity.
    mesh->libmesh_assert_valid_parallel_ids();

    // Make an adjacency data structure.
    DataTransferKit::LibmeshAdjacencies adjacencies( mesh );

    // Check the node adjacencies of the first hex elem.
    unsigned int num_nodes = 8;
    Teuchos::Array<Teuchos::Ptr<libMesh::Node>> elem_1_nodes;
    adjacencies.getLibmeshAdjacencies( Teuchos::ptr( hex_elem_1 ),
                                       elem_1_nodes );
    TEST_EQUALITY( num_nodes, elem_1_nodes.size() );
    for ( unsigned int n = 0; n < num_nodes; ++n )
    {
        TEST_EQUALITY( nodes[n]->id(), elem_1_nodes[n]->id() );
    };

    // Check the node adjacencies of the second hex elem.
    Teuchos::Array<Teuchos::Ptr<libMesh::Node>> elem_2_nodes;
    adjacencies.getLibmeshAdjacencies( Teuchos::ptr( hex_elem_2 ),
                                       elem_2_nodes );
    TEST_EQUALITY( num_nodes, elem_2_nodes.size() );
    for ( unsigned int n = 0; n < num_nodes; ++n )
    {
        TEST_EQUALITY( nodes[n]->id(), elem_2_nodes[n]->id() );
    };

    // Check the elem adjacencies of the first hex.
    Teuchos::Array<Teuchos::Ptr<libMesh::Elem>> elem_1_elems;
    adjacencies.getLibmeshAdjacencies( Teuchos::ptr( hex_elem_1 ),
                                       elem_1_elems );
    TEST_EQUALITY( 0, elem_1_elems.size() );

    // Check the elem adjacencies of the second hex.
    Teuchos::Array<Teuchos::Ptr<libMesh::Elem>> elem_2_elems;
    adjacencies.getLibmeshAdjacencies( Teuchos::ptr( hex_elem_2 ),
                                       elem_2_elems );
    TEST_EQUALITY( 0, elem_2_elems.size() );

    // Check the adjacencies of the nodes.
    for ( unsigned int n = 0; n < num_nodes; ++n )
    {
        Teuchos::Array<Teuchos::Ptr<libMesh::Elem>> node_elems;
        adjacencies.getLibmeshAdjacencies( Teuchos::ptr( nodes[n] ),
                                           node_elems );
        TEST_EQUALITY( 2, node_elems.size() );
        TEST_EQUALITY( hex_elem_2->id(), node_elems[0]->id() );
        TEST_EQUALITY( hex_elem_1->id(), node_elems[1]->id() );

        Teuchos::Array<Teuchos::Ptr<libMesh::Node>> node_nodes;
        adjacencies.getLibmeshAdjacencies( Teuchos::ptr( nodes[n] ),
                                           node_nodes );
        TEST_EQUALITY( 0, node_nodes.size() );
    }
}

//---------------------------------------------------------------------------//
// end tstLibmeshAdjacencies.cpp
//---------------------------------------------------------------------------//
