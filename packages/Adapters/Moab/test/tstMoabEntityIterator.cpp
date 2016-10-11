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
 * \file tstMoabEntityIterator.cpp
 * \author Stuart R. Slattery
 * \brief MoabEntityIterator unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_MoabEntityExtraData.hpp>
#include <DTK_MoabEntityIterator.hpp>
#include <DTK_MoabEntityIteratorRange.hpp>
#include <DTK_MoabEntityPredicates.hpp>
#include <DTK_MoabHelpers.hpp>
#include <DTK_MoabMeshSetIndexer.hpp>

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

#include <moab/Core.hpp>
#include <moab/Interface.hpp>
#include <moab/ParallelComm.hpp>

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
TEUCHOS_UNIT_TEST( MoabEntityIterator, hex_8_test )
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
    Teuchos::RCP<moab::Interface> moab_mesh = Teuchos::rcp( new moab::Core() );
    Teuchos::RCP<moab::ParallelComm> parallel_mesh = Teuchos::rcp(
        new moab::ParallelComm( moab_mesh.getRawPtr(), raw_comm ) );

    // Create the nodes.
    moab::ErrorCode error = moab::MB_SUCCESS;
    Teuchos::Array<moab::EntityHandle> nodes( 8 );
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
    error = moab_mesh->create_element( moab::MBHEX, nodes.getRawPtr(), 8,
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

    // Index the sets in the mesh.
    Teuchos::RCP<DataTransferKit::MoabMeshSetIndexer> set_indexer =
        Teuchos::rcp(
            new DataTransferKit::MoabMeshSetIndexer( parallel_mesh ) );

    // Make a list of hexes.
    unsigned num_hex = 2;
    std::vector<moab::EntityHandle> hex_entities( num_hex, hex_entity );

    // Make an iterator for the hex.
    std::function<bool( DataTransferKit::Entity )> all_pred =
        [=]( DataTransferKit::Entity e ) { return true; };
    Teuchos::RCP<DataTransferKit::MoabEntityIteratorRange> iterator_range =
        Teuchos::rcp( new DataTransferKit::MoabEntityIteratorRange() );
    iterator_range->d_moab_entities = hex_entities;
    DataTransferKit::EntityIterator entity_iterator =
        DataTransferKit::MoabEntityIterator(
            iterator_range, parallel_mesh.ptr(), set_indexer.ptr(), all_pred );

    // Test the entity iterator.
    TEST_EQUALITY( entity_iterator.size(), num_hex );
    TEST_ASSERT( entity_iterator == entity_iterator.begin() );
    TEST_ASSERT( entity_iterator != entity_iterator.end() );

    // Test the first entity under the iterator with a pointer dereference.
    DataTransferKit::EntityId hex_id = 90343;
    DataTransferKit::MoabHelpers::getGlobalIds( *parallel_mesh, &hex_entity, 1,
                                                &hex_id );
    TEST_EQUALITY( hex_id, entity_iterator->id() );
    TEST_EQUALITY( comm->getRank(), entity_iterator->ownerRank() );
    TEST_EQUALITY( space_dim, entity_iterator->topologicalDimension() );
    TEST_EQUALITY( space_dim, entity_iterator->physicalDimension() );

    int set_1_id = set_indexer->getIndexFromMeshSet( entity_set_1 );
    int set_2_id = set_indexer->getIndexFromMeshSet( entity_set_2 );
    TEST_ASSERT( entity_iterator->inBlock( set_1_id ) );
    TEST_ASSERT( !entity_iterator->inBlock( set_2_id ) );
    TEST_ASSERT( entity_iterator->onBoundary( set_1_id ) );
    TEST_ASSERT( !entity_iterator->onBoundary( set_2_id ) );

    Teuchos::RCP<DataTransferKit::EntityExtraData> extra_data_1 =
        entity_iterator->extraData();
    TEST_EQUALITY(
        hex_entity,
        Teuchos::rcp_dynamic_cast<DataTransferKit::MoabEntityExtraData>(
            extra_data_1 )
            ->d_moab_entity );

    Teuchos::Tuple<double, 6> hex_bounds_1;
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
    TEST_EQUALITY( hex_id, ( *entity_iterator ).id() );
    TEST_EQUALITY( comm->getRank(), ( *entity_iterator ).ownerRank() );
    TEST_EQUALITY( space_dim, ( *entity_iterator ).topologicalDimension() );
    TEST_EQUALITY( space_dim, ( *entity_iterator ).physicalDimension() );

    TEST_ASSERT( ( *entity_iterator ).inBlock( set_1_id ) );
    TEST_ASSERT( !( *entity_iterator ).inBlock( set_2_id ) );
    TEST_ASSERT( ( *entity_iterator ).onBoundary( set_1_id ) );
    TEST_ASSERT( !( *entity_iterator ).onBoundary( set_2_id ) );

    Teuchos::RCP<DataTransferKit::EntityExtraData> extra_data_2 =
        ( *entity_iterator ).extraData();
    TEST_EQUALITY(
        hex_entity,
        Teuchos::rcp_dynamic_cast<DataTransferKit::MoabEntityExtraData>(
            extra_data_2 )
            ->d_moab_entity );

    Teuchos::Tuple<double, 6> hex_bounds_2;
    ( *entity_iterator ).boundingBox( hex_bounds_2 );
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
    DataTransferKit::MoabMeshSetPredicate set_1_pred( entity_set_1,
                                                      set_indexer );
    DataTransferKit::EntityIterator set_1_iterator =
        DataTransferKit::MoabEntityIterator(
            iterator_range, parallel_mesh.ptr(), set_indexer.ptr(),
            set_1_pred.getFunction() );
    TEST_EQUALITY( set_1_iterator.size(), num_hex );

    // Make an iterator with a part 2 predicate.
    DataTransferKit::MoabMeshSetPredicate set_2_pred( entity_set_2,
                                                      set_indexer );
    DataTransferKit::EntityIterator set_2_iterator =
        DataTransferKit::MoabEntityIterator(
            iterator_range, parallel_mesh.ptr(), set_indexer.ptr(),
            set_2_pred.getFunction() );
    TEST_EQUALITY( set_2_iterator.size(), 0 );
}

//---------------------------------------------------------------------------//
// end tstMoabEntityIterator.cpp
//---------------------------------------------------------------------------//
