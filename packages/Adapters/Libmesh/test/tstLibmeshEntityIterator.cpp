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
 * \file tstLibmeshEntityIterator.cpp
 * \author Stuart R. Slattery
 * \brief LibmeshEntityIterator unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_BasicEntityPredicates.hpp>

#include <DTK_LibmeshAdjacencies.hpp>
#include <DTK_LibmeshEntity.hpp>
#include <DTK_LibmeshEntityExtraData.hpp>
#include <DTK_LibmeshEntityIterator.hpp>

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
TEUCHOS_UNIT_TEST( LibmeshEntityIterator, hex_8_test )
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

    // Make an adjacency data structure.
    DataTransferKit::LibmeshAdjacencies adjacencies( mesh );

    // Make an iterator for the hex.
    std::function<bool( DataTransferKit::Entity )> all_pred =
        [=]( DataTransferKit::Entity ) { return true; };
    DataTransferKit::EntityIterator entity_iterator =
        DataTransferKit::LibmeshEntityIterator<
            libMesh::Mesh::const_element_iterator>(
            mesh->elements_begin(), mesh->elements_begin(),
            mesh->elements_end(), mesh.ptr(),
            Teuchos::ptrFromRef( adjacencies ), all_pred );

    // Test the entity iterator.
    unsigned int num_hex = 1;
    TEST_EQUALITY( entity_iterator.size(), num_hex );
    TEST_ASSERT( entity_iterator == entity_iterator.begin() );
    TEST_ASSERT( entity_iterator != entity_iterator.end() );

    // Test the first entity under the iterator with a pointer dereference.
    TEST_EQUALITY( hex_elem->id(), entity_iterator->id() );
    TEST_EQUALITY( comm->getRank(), entity_iterator->ownerRank() );
    TEST_EQUALITY( space_dim, entity_iterator->topologicalDimension() );
    TEST_EQUALITY( space_dim, entity_iterator->physicalDimension() );

    TEST_ASSERT( entity_iterator->inBlock( subdomain_1_id ) );
    TEST_ASSERT( !entity_iterator->inBlock( subdomain_2_id ) );
    TEST_ASSERT( entity_iterator->onBoundary( boundary_1_id ) );
    TEST_ASSERT( !entity_iterator->onBoundary( boundary_2_id ) );

    Teuchos::RCP<DataTransferKit::EntityExtraData> extra_data_1 =
        entity_iterator->extraData();
    TEST_EQUALITY( hex_elem,
                   Teuchos::rcp_dynamic_cast<
                       DataTransferKit::LibmeshEntityExtraData<libMesh::Elem>>(
                       extra_data_1 )
                       ->d_libmesh_geom.getRawPtr() );

    Teuchos::Tuple<double, 6> hex_bounds_1;
    entity_iterator->boundingBox( hex_bounds_1 );
    TEST_EQUALITY( 0.0, hex_bounds_1[0] );
    TEST_EQUALITY( 0.0, hex_bounds_1[1] );
    TEST_EQUALITY( 0.0, hex_bounds_1[2] );
    TEST_EQUALITY( 1.0, hex_bounds_1[3] );
    TEST_EQUALITY( 1.0, hex_bounds_1[4] );
    TEST_EQUALITY( 1.0, hex_bounds_1[5] );

    // Test the end of the iterator.
    entity_iterator++;
    TEST_ASSERT( entity_iterator != entity_iterator.begin() );
    TEST_ASSERT( entity_iterator == entity_iterator.end() );

    // Make an iterator with a subdomain 1 predicate.
    DataTransferKit::BlockPredicate subdomain_1_pred(
        Teuchos::Array<int>( 1, subdomain_1_id ) );
    DataTransferKit::EntityIterator subdomain_1_iterator =
        DataTransferKit::LibmeshEntityIterator<
            libMesh::Mesh::const_element_iterator>(
            mesh->elements_begin(), mesh->elements_begin(),
            mesh->elements_end(), mesh.ptr(),
            Teuchos::ptrFromRef( adjacencies ),
            subdomain_1_pred.getFunction() );
    TEST_EQUALITY( subdomain_1_iterator.size(), num_hex );

    // Make an iterator with a subdomain 2 predicate.
    DataTransferKit::BlockPredicate subdomain_2_pred(
        Teuchos::Array<int>( 2, subdomain_2_id ) );
    DataTransferKit::EntityIterator subdomain_2_iterator =
        DataTransferKit::LibmeshEntityIterator<
            libMesh::Mesh::const_element_iterator>(
            mesh->elements_begin(), mesh->elements_begin(),
            mesh->elements_end(), mesh.ptr(),
            Teuchos::ptrFromRef( adjacencies ),
            subdomain_2_pred.getFunction() );
    TEST_EQUALITY( subdomain_2_iterator.size(), 0 );

    // Make a boundary iterator for the elems.
    DataTransferKit::BoundaryPredicate boundary_1_elem_pred(
        Teuchos::Array<int>( 1, subdomain_1_id ) );
    DataTransferKit::EntityIterator elem_boundary_it_1 =
        DataTransferKit::LibmeshEntityIterator<
            libMesh::Mesh::const_element_iterator>(
            mesh->elements_begin(), mesh->elements_begin(),
            mesh->elements_end(), mesh.ptr(),
            Teuchos::ptrFromRef( adjacencies ),
            boundary_1_elem_pred.getFunction() );
    TEST_EQUALITY( 1, elem_boundary_it_1.size() );

    // Make a boundary iterator for the elems.
    DataTransferKit::BoundaryPredicate boundary_2_elem_pred(
        Teuchos::Array<int>( 1, subdomain_2_id ) );
    DataTransferKit::EntityIterator elem_boundary_it_2 =
        DataTransferKit::LibmeshEntityIterator<
            libMesh::Mesh::const_element_iterator>(
            mesh->elements_begin(), mesh->elements_begin(),
            mesh->elements_end(), mesh.ptr(),
            Teuchos::ptrFromRef( adjacencies ),
            boundary_2_elem_pred.getFunction() );
    TEST_EQUALITY( 0, elem_boundary_it_2.size() );

    // Make an iterator for the nodes.
    DataTransferKit::EntityIterator node_iterator =
        DataTransferKit::LibmeshEntityIterator<
            libMesh::Mesh::const_node_iterator>(
            mesh->nodes_begin(), mesh->nodes_begin(), mesh->nodes_end(),
            mesh.ptr(), Teuchos::ptrFromRef( adjacencies ), all_pred );
    TEST_EQUALITY( 8, node_iterator.size() );

    // Make a boundary iterator for the nodes.
    DataTransferKit::BoundaryPredicate boundary_1_node_pred(
        Teuchos::Array<int>( 1, subdomain_1_id ) );
    DataTransferKit::EntityIterator node_boundary_it_1 =
        DataTransferKit::LibmeshEntityIterator<
            libMesh::Mesh::const_node_iterator>(
            mesh->nodes_begin(), mesh->nodes_begin(), mesh->nodes_end(),
            mesh.ptr(), Teuchos::ptrFromRef( adjacencies ),
            boundary_1_node_pred.getFunction() );
    TEST_EQUALITY( 0, node_boundary_it_1.size() );

    // Make a boundary iterator for the nodes.
    DataTransferKit::BoundaryPredicate boundary_2_node_pred(
        Teuchos::Array<int>( 1, subdomain_2_id ) );
    DataTransferKit::EntityIterator node_boundary_it_2 =
        DataTransferKit::LibmeshEntityIterator<
            libMesh::Mesh::const_node_iterator>(
            mesh->nodes_begin(), mesh->nodes_begin(), mesh->nodes_end(),
            mesh.ptr(), Teuchos::ptrFromRef( adjacencies ),
            boundary_2_node_pred.getFunction() );
    TEST_EQUALITY( 1, node_boundary_it_2.size() );
}

//---------------------------------------------------------------------------//
// end tstLibmeshEntityIterator.cpp
//---------------------------------------------------------------------------//
