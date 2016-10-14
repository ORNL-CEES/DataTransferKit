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

    // Make a range for the iterators.
    Teuchos::RCP<DataTransferKit::MoabEntityIteratorRange> iterator_range =
        Teuchos::rcp( new DataTransferKit::MoabEntityIteratorRange() );
    iterator_range->d_moab_entities = hex_entities;

    // Test the predicate for set 1.
    DataTransferKit::MoabMeshSetPredicate set_1_pred( entity_set_1,
                                                      set_indexer );
    DataTransferKit::EntityIterator set_1_iterator =
        DataTransferKit::MoabEntityIterator(
            iterator_range, parallel_mesh.ptr(), set_indexer.ptr(),
            set_1_pred.getFunction() );
    TEST_EQUALITY( set_1_iterator.size(), num_hex );

    // Test the predicate for set 2.
    DataTransferKit::MoabMeshSetPredicate set_2_pred( entity_set_2,
                                                      set_indexer );
    DataTransferKit::EntityIterator set_2_iterator =
        DataTransferKit::MoabEntityIterator(
            iterator_range, parallel_mesh.ptr(), set_indexer.ptr(),
            set_2_pred.getFunction() );
    TEST_EQUALITY( set_2_iterator.size(), 0 );

    // Test the vector predicate for part 1.
    std::vector<moab::EntityHandle> p1_vec( 1, entity_set_1 );
    DataTransferKit::MoabMeshSetPredicate entity_set_1_vec_pred( p1_vec,
                                                                 set_indexer );
    DataTransferKit::EntityIterator entity_set_1_vec_iterator =
        DataTransferKit::MoabEntityIterator(
            iterator_range, parallel_mesh.ptr(), set_indexer.ptr(),
            entity_set_1_vec_pred.getFunction() );
    TEST_EQUALITY( entity_set_1_vec_iterator.size(), num_hex );

    // Test the vector predicate for part 2.
    std::vector<moab::EntityHandle> p2_vec( 2, entity_set_2 );
    DataTransferKit::MoabMeshSetPredicate entity_set_2_vec_pred( p2_vec,
                                                                 set_indexer );
    DataTransferKit::EntityIterator entity_set_2_vec_iterator =
        DataTransferKit::MoabEntityIterator(
            iterator_range, parallel_mesh.ptr(), set_indexer.ptr(),
            entity_set_2_vec_pred.getFunction() );
    TEST_EQUALITY( entity_set_2_vec_iterator.size(), 0 );

    // Test a vector with 2 part 1's.
    std::vector<moab::EntityHandle> p11_vec( 2, entity_set_1 );
    DataTransferKit::MoabMeshSetPredicate entity_set_11_vec_pred( p11_vec,
                                                                  set_indexer );
    DataTransferKit::EntityIterator entity_set_11_vec_iterator =
        DataTransferKit::MoabEntityIterator(
            iterator_range, parallel_mesh.ptr(), set_indexer.ptr(),
            entity_set_11_vec_pred.getFunction() );
    TEST_EQUALITY( entity_set_11_vec_iterator.size(), num_hex );

    // Test a vector with a part 1 and part 2
    std::vector<moab::EntityHandle> p12_vec( 2 );
    p12_vec[0] = entity_set_1;
    p12_vec[1] = entity_set_2;
    DataTransferKit::MoabMeshSetPredicate entity_set_12_vec_pred( p12_vec,
                                                                  set_indexer );
    DataTransferKit::EntityIterator entity_set_12_vec_iterator =
        DataTransferKit::MoabEntityIterator(
            iterator_range, parallel_mesh.ptr(), set_indexer.ptr(),
            entity_set_12_vec_pred.getFunction() );
    TEST_EQUALITY( entity_set_12_vec_iterator.size(), 0 );
}

//---------------------------------------------------------------------------//
// end tstMoabEntityIterator.cpp
//---------------------------------------------------------------------------//
