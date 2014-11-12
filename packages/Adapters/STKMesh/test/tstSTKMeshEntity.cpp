//---------------------------------------------------------------------------//
/*!
 * \file tstSTKMeshEntity.cpp
 * \author Stuart R. Slattery
 * \brief STKMeshEntity unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_STKMeshEntity.hpp>
#include <DTK_STKMeshEntityExtraData.hpp>

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
TEUCHOS_UNIT_TEST( STKMeshEntity, hex_8_test )
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
    int part_1_id = part_1.id();
    std::string p2_name = "part_2";
    stk::mesh::Part& part_2 = meta_data.declare_part( p2_name );
    int part_2_id = part_2.id();

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
	bulk_data->declare_entity( stk::topology::ELEM_RANK, 23, part_1 );
    int num_nodes = 8;
    Teuchos::Array<stk::mesh::EntityId> node_ids( num_nodes );
    Teuchos::Array<stk::mesh::Entity> nodes( num_nodes );
    for ( int i = 0; i < num_nodes; ++i )
    {
	node_ids[i] = i + 5;
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

    // Create a DTK entity for the hex.
    DataTransferKit::Entity dtk_entity = 
	DataTransferKit::STKMeshEntity( hex_entity, bulk_data.ptr() );

    // Test the entity.
    TEST_EQUALITY( DataTransferKit::ENTITY_TYPE_VOLUME, dtk_entity.entityType() );
    TEST_EQUALITY( hex_id, dtk_entity.id() );
    TEST_EQUALITY( comm_rank, dtk_entity.ownerRank() );
    TEST_EQUALITY( space_dim, dtk_entity.physicalDimension() );
    TEST_ASSERT( dtk_entity.inBlock(part_1_id) );
    TEST_ASSERT( !dtk_entity.inBlock(part_2_id) );
    TEST_ASSERT( dtk_entity.onBoundary(part_1_id) );
    TEST_ASSERT( !dtk_entity.onBoundary(part_2_id) );

    Teuchos::RCP<DataTransferKit::EntityExtraData> extra_data =
	dtk_entity.extraData();
    TEST_EQUALITY( hex_entity.m_value,
		   Teuchos::rcp_dynamic_cast<DataTransferKit::STKMeshEntityExtraData>(
		       extra_data)->d_stk_entity.m_value );

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
// end tstSTKMeshEntity.cpp
//---------------------------------------------------------------------------//

