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
 * \file tstSTKMeshEntitySet.cpp
 * \author Stuart R. Slattery
 * \brief STKMeshEntitySet unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_STKMeshEntitySet.hpp>
#include <DTK_STKMeshEntityExtraData.hpp>
#include <DTK_EntitySet.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_topology/topology.hpp>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Hex-8 test.
TEUCHOS_UNIT_TEST( STKMeshEntitySet, hex_8_test )
{
    // Extract the raw mpi communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm =
        Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm<int> >( comm );
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm =
        mpi_comm->getRawMpiComm();
    MPI_Comm raw_comm = (*opaque_comm)();

    // Create meta data.
    int space_dim = 3;
    stk::mesh::MetaData meta_data( space_dim );

    // Make two parts.
    std::string p1_name = "part_1";
    stk::mesh::Part& part_1 = meta_data.declare_part( p1_name );
    stk::mesh::set_topology( part_1, stk::topology::HEX_8 );
    int part_1_id = part_1.mesh_meta_data_ordinal();
    std::string p2_name = "part_2";
    stk::mesh::Part& part_2 = meta_data.declare_part( p2_name );
    int part_2_id = part_2.mesh_meta_data_ordinal();

    // Make a coordinate field.
    stk::mesh::Field<double, stk::mesh::Cartesian3d>& coord_field =
        meta_data.declare_field<
        stk::mesh::Field<double, stk::mesh::Cartesian3d> >(
            stk::topology::NODE_RANK, "coordinates");
    meta_data.set_coordinate_field( &coord_field );
    stk::mesh::put_field( coord_field, part_1 );
    meta_data.commit();

    // Create bulk data.
    Teuchos::RCP<stk::mesh::BulkData> bulk_data =
        Teuchos::rcp( new stk::mesh::BulkData(meta_data,raw_comm) );
    bulk_data->modification_begin();

    // Make a hex-8.
    int comm_rank = comm->getRank();
    stk::mesh::EntityId hex_id = 23 + comm_rank;
    stk::mesh::Entity hex_entity =
        bulk_data->declare_entity( stk::topology::ELEM_RANK, hex_id, part_1 );
    unsigned num_nodes = 8;
    Teuchos::Array<stk::mesh::EntityId> node_ids( num_nodes );
    Teuchos::Array<stk::mesh::Entity> nodes( num_nodes );
    for ( unsigned i = 0; i < num_nodes; ++i )
    {
        node_ids[i] = num_nodes*comm_rank + i + 5;
        nodes[i] = bulk_data->declare_entity(
            stk::topology::NODE_RANK, node_ids[i], part_1 );
        bulk_data->declare_relation( hex_entity, nodes[i], i );
    }
    bulk_data->modification_end();

    // Create the node coordinates.
    double* node_coords = 0;
    node_coords = stk::mesh::field_data( coord_field, nodes[0] );
    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[1] );
    node_coords[0] = 1.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[2] );
    node_coords[0] = 1.0;
    node_coords[1] = 1.0;
    node_coords[2] = 0.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[3] );
    node_coords[0] = 0.0;
    node_coords[1] = 1.0;
    node_coords[2] = 0.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[4] );
    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = 1.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[5] );
    node_coords[0] = 1.0;
    node_coords[1] = 0.0;
    node_coords[2] = 1.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[6] );
    node_coords[0] = 1.0;
    node_coords[1] = 1.0;
    node_coords[2] = 1.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[7] );
    node_coords[0] = 0.0;
    node_coords[1] = 1.0;
    node_coords[2] = 1.0;

    // Create an entity set.
    Teuchos::RCP<DataTransferKit::EntitySet> entity_set =
        Teuchos::rcp( new DataTransferKit::STKMeshEntitySet(bulk_data) );

    // Test the set.
    Teuchos::RCP<const Teuchos::Comm<int> > set_comm =
        entity_set->communicator();
    TEST_EQUALITY( set_comm->getRank(), comm->getRank() );
    TEST_EQUALITY( set_comm->getSize(), comm->getSize() );
    TEST_EQUALITY( space_dim, entity_set->physicalDimension() );

    // Make an iterator for the hex.
    std::function<bool(DataTransferKit::Entity)> all_pred =
        [=] (DataTransferKit::Entity e){return true;};
    DataTransferKit::EntityIterator volume_iterator =
        entity_set->entityIterator( space_dim, all_pred );

    // Test the volume iterator.
    TEST_EQUALITY( volume_iterator.size(), 1 );
    TEST_ASSERT( volume_iterator == volume_iterator.begin() );
    TEST_ASSERT( volume_iterator != volume_iterator.end() );

    // Test the volume under the iterator.
    TEST_EQUALITY( hex_id, volume_iterator->id() );
    TEST_EQUALITY( comm_rank, volume_iterator->ownerRank() );
    TEST_EQUALITY( space_dim, volume_iterator->topologicalDimension() );
    TEST_EQUALITY( space_dim, volume_iterator->physicalDimension() );
    TEST_ASSERT( volume_iterator->inBlock(part_1_id) );
    TEST_ASSERT( !volume_iterator->inBlock(part_2_id) );
    TEST_ASSERT( volume_iterator->onBoundary(part_1_id) );
    TEST_ASSERT( !volume_iterator->onBoundary(part_2_id) );

    Teuchos::RCP<DataTransferKit::EntityExtraData> extra_data_1 =
        volume_iterator->extraData();
    TEST_EQUALITY( hex_entity,
                   Teuchos::rcp_dynamic_cast<DataTransferKit::STKMeshEntityExtraData>(
                       extra_data_1)->d_stk_entity );

    Teuchos::Tuple<double,6> hex_bounds_1;
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
    TEST_EQUALITY( node_iterator.size(), num_nodes );
    TEST_ASSERT( node_iterator == node_iterator.begin() );
    TEST_ASSERT( node_iterator != node_iterator.end() );
    DataTransferKit::EntityIterator node_begin = node_iterator.begin();
    DataTransferKit::EntityIterator node_end = node_iterator.end();
    auto node_id_it = node_ids.begin();
    for ( node_iterator = node_begin;
          node_iterator != node_end;
          ++node_iterator, ++node_id_it )
    {
        TEST_EQUALITY( node_iterator->id(), *node_id_it );
    }

    // Get each entity and check.
    DataTransferKit::Entity set_hex;
    entity_set->getEntity( hex_id, space_dim, set_hex );
    TEST_EQUALITY( set_hex.id(), hex_id );
    for ( unsigned i = 0; i < num_nodes; ++i )
    {
        DataTransferKit::Entity set_node;
        entity_set->getEntity( node_ids[i], 0, set_node );
        TEST_EQUALITY( set_node.id(), node_ids[i] );
    }

    // Check the adjacency function.
    Teuchos::Array<DataTransferKit::Entity> hex_adjacent_volumes;
    entity_set->getAdjacentEntities( set_hex,
                                     space_dim,
                                     hex_adjacent_volumes );
    TEST_EQUALITY( 0, hex_adjacent_volumes.size() );

    Teuchos::Array<DataTransferKit::Entity> hex_adjacent_nodes;
    entity_set->getAdjacentEntities( set_hex, 0, hex_adjacent_nodes );
    TEST_EQUALITY( num_nodes, hex_adjacent_nodes.size() );
    for ( unsigned i = 0; i < num_nodes; ++i )
    {
        TEST_EQUALITY( hex_adjacent_nodes[i].id(), node_ids[i] );
    }

    for ( unsigned i = 0; i < num_nodes; ++i )
    {
        Teuchos::Array<DataTransferKit::Entity> node_adjacent_volumes;
        entity_set->getAdjacentEntities( hex_adjacent_nodes[i],
                                         space_dim,
                                         node_adjacent_volumes );
        TEST_EQUALITY( 1, node_adjacent_volumes.size() );
        TEST_EQUALITY( node_adjacent_volumes[0].id(), hex_id );
    }
}

//---------------------------------------------------------------------------//
// end tstSTKMeshEntitySet.cpp
//---------------------------------------------------------------------------//
