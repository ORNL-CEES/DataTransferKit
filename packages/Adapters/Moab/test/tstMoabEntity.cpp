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

    // Test the entity.
    TEST_EQUALITY( DataTransferKit::ENTITY_TYPE_VOLUME, dtk_entity.entityType() );
    TEST_EQUALITY( hex_entity, dtk_entity.id() );
    TEST_EQUALITY( comm->getRank(), dtk_entity.ownerRank() );
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

