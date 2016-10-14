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
 * \file tstLibmeshEntitySet.cpp
 * \author Stuart R. Slattery
 * \brief LibmeshEntitySet unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_LibmeshEntityExtraData.hpp>
#include <DTK_LibmeshEntitySet.hpp>

#include <DTK_EntitySet.hpp>

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
TEUCHOS_UNIT_TEST( LibmeshEntitySet, hex_8_test )
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

    // Create an entity set.
    Teuchos::RCP<DataTransferKit::EntitySet> entity_set =
        Teuchos::rcp( new DataTransferKit::LibmeshEntitySet( mesh ) );

    // Test the set.
    Teuchos::RCP<const Teuchos::Comm<int>> set_comm =
        entity_set->communicator();
    TEST_EQUALITY( set_comm->getRank(), comm->getRank() );
    TEST_EQUALITY( set_comm->getSize(), comm->getSize() );
    TEST_EQUALITY( space_dim, entity_set->physicalDimension() );

    // Make an iterator for the hex.
    std::function<bool( DataTransferKit::Entity )> all_pred =
        [=]( DataTransferKit::Entity ) { return true; };
    DataTransferKit::EntityIterator volume_iterator =
        entity_set->entityIterator( 3, all_pred );

    // Test the volume iterator.
    TEST_EQUALITY( volume_iterator.size(), 1 );
    TEST_ASSERT( volume_iterator == volume_iterator.begin() );
    TEST_ASSERT( volume_iterator != volume_iterator.end() );

    // Test the volume under the iterator.
    TEST_EQUALITY( hex_elem->id(), volume_iterator->id() );
    TEST_EQUALITY( comm->getRank(), volume_iterator->ownerRank() );
    TEST_EQUALITY( space_dim, volume_iterator->topologicalDimension() );
    TEST_EQUALITY( space_dim, volume_iterator->physicalDimension() );

    Teuchos::RCP<DataTransferKit::EntityExtraData> extra_data_1 =
        volume_iterator->extraData();
    TEST_EQUALITY( hex_elem,
                   Teuchos::rcp_dynamic_cast<
                       DataTransferKit::LibmeshEntityExtraData<libMesh::Elem>>(
                       extra_data_1 )
                       ->d_libmesh_geom.getRawPtr() );

    Teuchos::Tuple<double, 6> hex_bounds_1;
    volume_iterator->boundingBox( hex_bounds_1 );
    TEST_EQUALITY( 0.0, hex_bounds_1[0] );
    TEST_EQUALITY( 0.0, hex_bounds_1[1] );
    TEST_EQUALITY( 0.0, hex_bounds_1[2] );
    TEST_EQUALITY( 1.0, hex_bounds_1[3] );
    TEST_EQUALITY( 1.0, hex_bounds_1[4] );
    TEST_EQUALITY( 1.0, hex_bounds_1[5] );

    // Test the end of the iterator.
    volume_iterator++;
    TEST_ASSERT( volume_iterator != volume_iterator.begin() );
    TEST_ASSERT( volume_iterator == volume_iterator.end() );

    // Make an iterator for the nodes.
    DataTransferKit::EntityIterator node_iterator =
        entity_set->entityIterator( 0, all_pred );

    // Test the node iterator.
    unsigned num_nodes = 8;
    TEST_EQUALITY( node_iterator.size(), num_nodes );
    TEST_ASSERT( node_iterator == node_iterator.begin() );
    TEST_ASSERT( node_iterator != node_iterator.end() );
    DataTransferKit::EntityIterator node_begin = node_iterator.begin();
    DataTransferKit::EntityIterator node_end = node_iterator.end();
    auto node_id_it = nodes.begin();
    for ( node_iterator = node_begin; node_iterator != node_end;
          ++node_iterator, ++node_id_it )
    {
        TEST_EQUALITY( node_iterator->id(), ( *node_id_it )->id() );
    }

    // Get each entity and check.
    DataTransferKit::Entity set_hex;
    entity_set->getEntity( hex_elem->id(), 3, set_hex );
    TEST_EQUALITY( set_hex.id(), hex_elem->id() );
    for ( unsigned i = 0; i < num_nodes; ++i )
    {
        DataTransferKit::Entity set_node;
        entity_set->getEntity( nodes[i]->id(), 0, set_node );
        TEST_EQUALITY( set_node.id(), nodes[i]->id() );
    }

    // Check the adjacency function.
    Teuchos::Array<DataTransferKit::Entity> hex_adjacent_volumes;
    entity_set->getAdjacentEntities( set_hex, 3, hex_adjacent_volumes );
    TEST_EQUALITY( 0, hex_adjacent_volumes.size() );

    Teuchos::Array<DataTransferKit::Entity> hex_adjacent_nodes;
    entity_set->getAdjacentEntities( set_hex, 0, hex_adjacent_nodes );
    TEST_EQUALITY( num_nodes, hex_adjacent_nodes.size() );
    for ( unsigned i = 0; i < num_nodes; ++i )
    {
        TEST_EQUALITY( hex_adjacent_nodes[i].id(), nodes[i]->id() );
    }

    for ( unsigned i = 0; i < num_nodes; ++i )
    {
        Teuchos::Array<DataTransferKit::Entity> node_adjacent_volumes;
        entity_set->getAdjacentEntities( hex_adjacent_nodes[i], 3,
                                         node_adjacent_volumes );
        TEST_EQUALITY( 1, node_adjacent_volumes.size() );
        TEST_EQUALITY( node_adjacent_volumes[0].id(), hex_elem->id() );
    }
}

//---------------------------------------------------------------------------//
// end tstLibmeshEntitySet.cpp
//---------------------------------------------------------------------------//
