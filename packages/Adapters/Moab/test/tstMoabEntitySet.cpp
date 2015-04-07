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
 * \file tstMoabEntitySet.cpp
 * \author Stuart R. Slattery
 * \brief MoabEntitySet unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_MoabEntitySet.hpp>
#include <DTK_MoabEntityExtraData.hpp>
#include <DTK_MoabMeshSetIndexer.hpp>
#include <DTK_MoabHelpers.hpp>
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

#include <MBInterface.hpp>
#include <MBParallelComm.hpp>
#include <MBCore.hpp>

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
TEUCHOS_UNIT_TEST( MoabEntitySet, hex_8_test )
{
    // Extract the raw mpi communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm = 
	Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm<int> >( comm );
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm = 
	mpi_comm->getRawMpiComm();
    MPI_Comm raw_comm = (*opaque_comm)();

    // Create the mesh.
    int space_dim = 3;
    Teuchos::RCP<moab::Interface> moab_mesh = Teuchos::rcp( new moab::Core() );
    Teuchos::RCP<moab::ParallelComm> parallel_mesh =
	Teuchos::rcp( new moab::ParallelComm(moab_mesh.getRawPtr(),raw_comm) );

    // Create the nodes.
    moab::ErrorCode error = moab::MB_SUCCESS;
    Teuchos::Array<moab::EntityHandle> nodes(8);
    double node_coords[3];
    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;
    error = moab_mesh->create_vertex( node_coords, nodes[0] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 1.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;
    error = moab_mesh->create_vertex( node_coords, nodes[1] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 1.0;
    node_coords[1] = 1.0;
    node_coords[2] = 0.0;
    error = moab_mesh->create_vertex( node_coords, nodes[2] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 0.0;
    node_coords[1] = 1.0;
    node_coords[2] = 0.0;
    error = moab_mesh->create_vertex( node_coords, nodes[3] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = 1.0;
    error = moab_mesh->create_vertex( node_coords, nodes[4] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 1.0;
    node_coords[1] = 0.0;
    node_coords[2] = 1.0;
    error = moab_mesh->create_vertex( node_coords, nodes[5] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 1.0;
    node_coords[1] = 1.0;
    node_coords[2] = 1.0;
    error = moab_mesh->create_vertex( node_coords, nodes[6] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 0.0;
    node_coords[1] = 1.0;
    node_coords[2] = 1.0;
    error = moab_mesh->create_vertex( node_coords, nodes[7] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    // Make a hex-8.
    moab::EntityHandle hex_entity;
    error = moab_mesh->create_element( moab::MBHEX,
				       nodes.getRawPtr(),
				       8,
				       hex_entity );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    // Make 2 entity sets.
    moab::EntityHandle entity_set_1;
    error = moab_mesh->create_meshset( 0, entity_set_1 );
    TEST_EQUALITY( error, moab::MB_SUCCESS );
    moab::EntityHandle entity_set_2;
    error = moab_mesh->create_meshset( 0, entity_set_2 );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    // Put the hex-8 in the first entity set.
    error = moab_mesh->add_entities( entity_set_1, &hex_entity, 1 );

    // Create an entity set.
    Teuchos::RCP<DataTransferKit::MoabMeshSetIndexer> set_indexer = Teuchos::rcp(
	new DataTransferKit::MoabMeshSetIndexer(parallel_mesh) );
    Teuchos::RCP<DataTransferKit::EntitySet> entity_set =
	Teuchos::rcp( new DataTransferKit::MoabEntitySet(
			  parallel_mesh,set_indexer) );

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
	entity_set->entityIterator( 3, all_pred );

    // Test the volume iterator.
    TEST_EQUALITY( volume_iterator.size(), 1 );
    TEST_ASSERT( volume_iterator == volume_iterator.begin() );
    TEST_ASSERT( volume_iterator != volume_iterator.end() );

    // Test the volume under the iterator.
    DataTransferKit::EntityId hex_id = 90343;
    DataTransferKit::MoabHelpers::getGlobalIds(
	*parallel_mesh, &hex_entity, 1, &hex_id );
    TEST_EQUALITY( hex_id, volume_iterator->id() );
    TEST_EQUALITY( comm->getRank(), volume_iterator->ownerRank() );
    TEST_EQUALITY( space_dim, volume_iterator->topologicalDimension() );
    TEST_EQUALITY( space_dim, volume_iterator->physicalDimension() );

    Teuchos::RCP<DataTransferKit::EntityExtraData> extra_data_1 =
	volume_iterator->extraData();
    TEST_EQUALITY( hex_entity,
		   Teuchos::rcp_dynamic_cast<DataTransferKit::MoabEntityExtraData>(
		       extra_data_1)->d_moab_entity );

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
    unsigned num_nodes = 8;
    std::vector<DataTransferKit::EntityId> node_ids( num_nodes );
    DataTransferKit::MoabHelpers::getGlobalIds( *parallel_mesh,
						nodes.getRawPtr(),
						num_nodes,
						node_ids.data() );
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
    entity_set->getEntity( hex_id, 3, set_hex );
    TEST_EQUALITY( set_hex.id(), hex_id );
    for ( unsigned i = 0; i < num_nodes; ++i )
    {
	DataTransferKit::Entity set_node;
	entity_set->getEntity( node_ids[i], 0, set_node );
	TEST_EQUALITY( set_node.id(), node_ids[i] );
    }

    // Check the adjacency function.
    Teuchos::Array<DataTransferKit::Entity> hex_adjacent_volumes;
    entity_set->getAdjacentEntities( set_hex, 3,
				     hex_adjacent_volumes );
    TEST_EQUALITY( 0, hex_adjacent_volumes.size() );

    Teuchos::Array<DataTransferKit::Entity> hex_adjacent_nodes;
    entity_set->getAdjacentEntities( set_hex, 0,
				     hex_adjacent_nodes );
    TEST_EQUALITY( num_nodes, hex_adjacent_nodes.size() );
    for ( unsigned i = 0; i < num_nodes; ++i )
    {
	TEST_EQUALITY( hex_adjacent_nodes[i].id(), node_ids[i] );
    }

    for ( unsigned i = 0; i < num_nodes; ++i )
    {
	Teuchos::Array<DataTransferKit::Entity> node_adjacent_volumes;
	entity_set->getAdjacentEntities( hex_adjacent_nodes[i], 3,
					 node_adjacent_volumes );
	TEST_EQUALITY( 1, node_adjacent_volumes.size() );
	TEST_EQUALITY( node_adjacent_volumes[0].id(), hex_id );
    }
}

//---------------------------------------------------------------------------//
// end tstMoabEntitySet.cpp
//---------------------------------------------------------------------------//
