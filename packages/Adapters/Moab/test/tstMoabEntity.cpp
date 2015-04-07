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
 * \file tstMoabEntity.cpp
 * \author Stuart R. Slattery
 * \brief MoabEntity unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_MoabEntity.hpp>
#include <DTK_MoabEntityExtraData.hpp>
#include <DTK_MoabMeshSetIndexer.hpp>
#include <DTK_MoabHelpers.hpp>

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
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>

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
TEUCHOS_UNIT_TEST( MoabEntity, hex_8_test )
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

    // Index the sets in the mesh.
    Teuchos::RCP<DataTransferKit::MoabMeshSetIndexer> set_indexer =
	Teuchos::rcp( new DataTransferKit::MoabMeshSetIndexer(parallel_mesh) );

    // Create a DTK entity for the hex.
    DataTransferKit::Entity dtk_entity = DataTransferKit::MoabEntity( 
	hex_entity, parallel_mesh.ptr(), set_indexer.ptr() );

    // Print out the entity.
    Teuchos::RCP<Teuchos::FancyOStream>
	fancy_out = Teuchos::VerboseObjectBase::getDefaultOStream();
    dtk_entity.describe( *fancy_out );

    // Test the entity.
    DataTransferKit::EntityId hex_id = 90343;
    DataTransferKit::MoabHelpers::getGlobalIds(
	*parallel_mesh, &hex_entity, 1, &hex_id );
    TEST_EQUALITY( hex_id, dtk_entity.id() );
    TEST_EQUALITY( comm->getRank(), dtk_entity.ownerRank() );
    TEST_EQUALITY( space_dim, dtk_entity.topologicalDimension() );
    TEST_EQUALITY( space_dim, dtk_entity.physicalDimension() );

    int set_1_id = set_indexer->getIndexFromMeshSet( entity_set_1 );
    int set_2_id = set_indexer->getIndexFromMeshSet( entity_set_2 );
    TEST_ASSERT( dtk_entity.inBlock(set_1_id) );
    TEST_ASSERT( !dtk_entity.inBlock(set_2_id) );
    TEST_ASSERT( dtk_entity.onBoundary(set_1_id) );
    TEST_ASSERT( !dtk_entity.onBoundary(set_2_id) );

    Teuchos::RCP<DataTransferKit::EntityExtraData> extra_data =
	dtk_entity.extraData();
    TEST_EQUALITY( hex_entity,
		   Teuchos::rcp_dynamic_cast<DataTransferKit::MoabEntityExtraData>(
		       extra_data)->d_moab_entity );

    Teuchos::Tuple<double,6> hex_bounds;
    dtk_entity.boundingBox( hex_bounds );
    TEST_EQUALITY( 0.0, hex_bounds[0] );
    TEST_EQUALITY( 0.0, hex_bounds[1] );
    TEST_EQUALITY( 0.0, hex_bounds[2] );
    TEST_EQUALITY( 1.0, hex_bounds[3] );
    TEST_EQUALITY( 1.0, hex_bounds[4] );
    TEST_EQUALITY( 1.0, hex_bounds[5] );
}

//---------------------------------------------------------------------------//
// end tstMoabEntity.cpp
//---------------------------------------------------------------------------//

