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
 * \file tstLibmeshEntity.cpp
 * \author Stuart R. Slattery
 * \brief LibmeshEntity unit tests.
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
#include <DTK_LibmeshEntity.hpp>
#include <DTK_LibmeshEntityExtraData.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>

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
    libMesh::Elem *hex_elem = mesh->add_elem( new libMesh::Hex8 );
    hex_elem->processor_id() = rank;
    hex_elem->set_id() = 2 * rank;
    for ( int i = 0; i < 8; ++i )
        hex_elem->set_node( i ) = nodes[i];

    // Make 2 subdomains and put the hex-8 in the first subdomain.
    int subdomain_1_id = 1;
    int subdomain_2_id = 2;
    std::set<libMesh::subdomain_id_type> subdomain_ids;
    subdomain_ids.insert( subdomain_1_id );
    subdomain_ids.insert( subdomain_2_id );
    hex_elem->subdomain_id() = subdomain_1_id;

    // Make 2 boundaries and add the first elem side to one and first node to
    // the second.
    int boundary_1_id = 1;
    int boundary_2_id = 2;
    mesh->get_boundary_info().add_side( hex_elem, 0, boundary_1_id );
    mesh->get_boundary_info().add_node( nodes[0], boundary_2_id );

    // Check libmesh validity.
    mesh->libmesh_assert_valid_parallel_ids();

    // Make an adjacency data structure.
    DataTransferKit::LibmeshAdjacencies adjacencies( mesh );

    // Create a DTK entity for the hex.
    DataTransferKit::Entity dtk_entity =
        DataTransferKit::LibmeshEntity<libMesh::Elem>(
            Teuchos::ptr( hex_elem ), mesh.ptr(),
            Teuchos::ptrFromRef( adjacencies ) );

    // Print out the entity.
    Teuchos::RCP<Teuchos::FancyOStream> fancy_out =
        Teuchos::VerboseObjectBase::getDefaultOStream();
    dtk_entity.describe( *fancy_out );

    // Test the entity.
    TEST_EQUALITY( hex_elem->id(), dtk_entity.id() );
    TEST_EQUALITY( rank, dtk_entity.ownerRank() );
    TEST_EQUALITY( space_dim, dtk_entity.topologicalDimension() );
    TEST_EQUALITY( space_dim, dtk_entity.physicalDimension() );

    TEST_ASSERT( dtk_entity.inBlock( subdomain_1_id ) );
    TEST_ASSERT( !dtk_entity.inBlock( subdomain_2_id ) );

    TEST_ASSERT( dtk_entity.onBoundary( boundary_1_id ) );
    TEST_ASSERT( !dtk_entity.onBoundary( boundary_2_id ) );

    Teuchos::RCP<DataTransferKit::EntityExtraData> elem_extra_data =
        dtk_entity.extraData();
    TEST_EQUALITY( hex_elem,
                   Teuchos::rcp_dynamic_cast<
                       DataTransferKit::LibmeshEntityExtraData<libMesh::Elem>>(
                       elem_extra_data )
                       ->d_libmesh_geom.getRawPtr() );

    Teuchos::Tuple<double, 6> hex_bounds;
    dtk_entity.boundingBox( hex_bounds );
    TEST_EQUALITY( 0.0, hex_bounds[0] );
    TEST_EQUALITY( 0.0, hex_bounds[1] );
    TEST_EQUALITY( 0.0, hex_bounds[2] );
    TEST_EQUALITY( 1.0, hex_bounds[3] );
    TEST_EQUALITY( 1.0, hex_bounds[4] );
    TEST_EQUALITY( 1.0, hex_bounds[5] );

    // Test a node.
    DataTransferKit::Entity dtk_node =
        DataTransferKit::LibmeshEntity<libMesh::Node>(
            Teuchos::ptr( nodes[0] ), mesh.ptr(),
            Teuchos::ptrFromRef( adjacencies ) );
    TEST_EQUALITY( nodes[0]->id(), dtk_node.id() );
    TEST_EQUALITY( rank, dtk_node.ownerRank() );
    TEST_EQUALITY( 0, dtk_node.topologicalDimension() );
    TEST_EQUALITY( space_dim, dtk_node.physicalDimension() );

    TEST_ASSERT( dtk_node.inBlock( subdomain_1_id ) );
    TEST_ASSERT( !dtk_node.inBlock( subdomain_2_id ) );

    TEST_ASSERT( !dtk_node.onBoundary( boundary_1_id ) );
    TEST_ASSERT( dtk_node.onBoundary( boundary_2_id ) );

    Teuchos::RCP<DataTransferKit::EntityExtraData> node_extra_data =
        dtk_node.extraData();
    TEST_EQUALITY( nodes[0],
                   Teuchos::rcp_dynamic_cast<
                       DataTransferKit::LibmeshEntityExtraData<libMesh::Node>>(
                       node_extra_data )
                       ->d_libmesh_geom.getRawPtr() );

    Teuchos::Tuple<double, 6> node_bounds;
    dtk_node.boundingBox( node_bounds );
    TEST_EQUALITY( 0.0, node_bounds[0] );
    TEST_EQUALITY( 0.0, node_bounds[1] );
    TEST_EQUALITY( 0.0, node_bounds[2] );
    TEST_EQUALITY( 0.0, node_bounds[3] );
    TEST_EQUALITY( 0.0, node_bounds[4] );
    TEST_EQUALITY( 0.0, node_bounds[5] );

    // Print out the node.
    dtk_node.describe( *fancy_out );
}

//---------------------------------------------------------------------------//
// end tstLibmeshEntity.cpp
//---------------------------------------------------------------------------//
