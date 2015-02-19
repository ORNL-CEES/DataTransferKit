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
 * \file tstSTKMeshEntityIterator.cpp
 * \author Stuart R. Slattery
 * \brief STKMeshEntityIterator unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_STKMeshEntityIterator.hpp>
#include <DTK_STKMeshEntityIteratorRange.hpp>
#include <DTK_STKMeshEntityExtraData.hpp>
#include <DTK_STKMeshEntityPredicates.hpp>

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
#include <stk_mesh/base/Selector.hpp>
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
TEUCHOS_UNIT_TEST( STKMeshEntityIterator, hex_8_test )
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
    int num_nodes = 8;
    Teuchos::Array<stk::mesh::EntityId> node_ids( num_nodes );
    Teuchos::Array<stk::mesh::Entity> nodes( num_nodes );
    for ( int i = 0; i < num_nodes; ++i )
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

    // Make a list of hexes.
    unsigned num_hex = 2;
    std::vector<stk::mesh::Entity> hex_entities( num_hex, hex_entity );
    
    // Make an iterator for the hex.
    std::function<bool(DataTransferKit::Entity)> all_pred = 
	[=] (DataTransferKit::Entity e){return true;};
    Teuchos::RCP<DataTransferKit::STKMeshEntityIteratorRange> iterator_range =
	Teuchos::rcp( new DataTransferKit::STKMeshEntityIteratorRange() );
    iterator_range->d_stk_entities = hex_entities;
    DataTransferKit::EntityIterator entity_iterator = 
	DataTransferKit::STKMeshEntityIterator(
	    iterator_range, bulk_data.ptr(), all_pred );

    // Test the entity iterator.
    TEST_EQUALITY( entity_iterator.size(), num_hex );
    TEST_ASSERT( entity_iterator == entity_iterator.begin() );
    TEST_ASSERT( entity_iterator != entity_iterator.end() );

    // Test the first entity under the iterator with a pointer dereference.
    TEST_EQUALITY( DataTransferKit::ENTITY_TYPE_VOLUME, entity_iterator->entityType() );
    TEST_EQUALITY( hex_id, entity_iterator->id() );
    TEST_EQUALITY( comm_rank, entity_iterator->ownerRank() );
    TEST_EQUALITY( space_dim, entity_iterator->physicalDimension() );
    TEST_ASSERT( entity_iterator->inBlock(part_1_id) );
    TEST_ASSERT( !entity_iterator->inBlock(part_2_id) );
    TEST_ASSERT( entity_iterator->onBoundary(part_1_id) );
    TEST_ASSERT( !entity_iterator->onBoundary(part_2_id) );

    Teuchos::RCP<DataTransferKit::EntityExtraData> extra_data_1 =
	entity_iterator->extraData();
    TEST_EQUALITY( hex_entity,
		   Teuchos::rcp_dynamic_cast<DataTransferKit::STKMeshEntityExtraData>(
		       extra_data_1)->d_stk_entity );

    Teuchos::Tuple<double,6> hex_bounds_1;
    entity_iterator->boundingBox( hex_bounds_1 );
    TEST_EQUALITY( 0.0, hex_bounds_1[0] );
    TEST_EQUALITY( 0.0, hex_bounds_1[1] );
    TEST_EQUALITY( 0.0, hex_bounds_1[2] );
    TEST_EQUALITY( 1.0, hex_bounds_1[3] );
    TEST_EQUALITY( 1.0, hex_bounds_1[4] );
    TEST_EQUALITY( 1.0, hex_bounds_1[5] );

    // Increment the iterator
    ++entity_iterator;

    // Test the second entity under the iterator with a reference dereference.
    TEST_EQUALITY( DataTransferKit::ENTITY_TYPE_VOLUME, 
		   (*entity_iterator).entityType() );
    TEST_EQUALITY( hex_id, (*entity_iterator).id() );
    TEST_EQUALITY( comm_rank, (*entity_iterator).ownerRank() );
    TEST_EQUALITY( space_dim, (*entity_iterator).physicalDimension() );
    TEST_ASSERT( (*entity_iterator).inBlock(part_1_id) );
    TEST_ASSERT( !(*entity_iterator).inBlock(part_2_id) );
    TEST_ASSERT( (*entity_iterator).onBoundary(part_1_id) );
    TEST_ASSERT( !(*entity_iterator).onBoundary(part_2_id) );

    Teuchos::RCP<DataTransferKit::EntityExtraData> extra_data_2 =
	(*entity_iterator).extraData();
    TEST_EQUALITY( hex_entity,
		   Teuchos::rcp_dynamic_cast<DataTransferKit::STKMeshEntityExtraData>(
		       extra_data_2)->d_stk_entity );

    Teuchos::Tuple<double,6> hex_bounds_2;
    (*entity_iterator).boundingBox( hex_bounds_2 );
    TEST_EQUALITY( 0.0, hex_bounds_2[0] );
    TEST_EQUALITY( 0.0, hex_bounds_2[1] );
    TEST_EQUALITY( 0.0, hex_bounds_2[2] );
    TEST_EQUALITY( 1.0, hex_bounds_2[3] );
    TEST_EQUALITY( 1.0, hex_bounds_2[4] );
    TEST_EQUALITY( 1.0, hex_bounds_2[5] );

    // Test the end of the iterator.
    entity_iterator++;
    TEST_ASSERT( entity_iterator != entity_iterator.begin() );
    TEST_ASSERT( entity_iterator == entity_iterator.end() );

    // Make an iterator with a part 1 predicate.
    stk::mesh::Selector select_1( part_1 );
    DataTransferKit::STKSelectorPredicate part_1_pred( select_1 );
    DataTransferKit::EntityIterator part_1_iterator =
	DataTransferKit::STKMeshEntityIterator(
	    iterator_range, bulk_data.ptr(), part_1_pred.getFunction() );
    TEST_EQUALITY( part_1_iterator.size(), num_hex );

    // Make an iterator with a part 2 predicate.
    stk::mesh::Selector select_2( part_2 );
    DataTransferKit::STKSelectorPredicate part_2_pred( select_2 );
    DataTransferKit::EntityIterator part_2_iterator =
	DataTransferKit::STKMeshEntityIterator(
	    iterator_range, bulk_data.ptr(), part_2_pred.getFunction() );
    TEST_EQUALITY( part_2_iterator.size(), 0 );
}

//---------------------------------------------------------------------------//
// end tstSTKMeshEntityIterator.cpp
//---------------------------------------------------------------------------//
