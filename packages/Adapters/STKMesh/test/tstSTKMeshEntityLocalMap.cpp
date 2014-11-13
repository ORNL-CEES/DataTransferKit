//---------------------------------------------------------------------------//
/*!
 * \file tstSTKMeshEntityLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief STKMeshEntityLocalMap unit tests.
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
#include <DTK_STKMeshEntityLocalMap.hpp>

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
    node_coords[0] = 2.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[2] );
    node_coords[0] = 2.0;
    node_coords[1] = 2.0;
    node_coords[2] = 0.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[3] );
    node_coords[0] = 0.0;
    node_coords[1] = 2.0;
    node_coords[2] = 0.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[4] );
    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = 2.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[5] );
    node_coords[0] = 2.0;
    node_coords[1] = 0.0;
    node_coords[2] = 2.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[6] );
    node_coords[0] = 2.0;
    node_coords[1] = 2.0;
    node_coords[2] = 2.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[7] );
    node_coords[0] = 0.0;
    node_coords[1] = 2.0;
    node_coords[2] = 2.0;
    
    // Create a local map from the bulk data.
    Teuchos::RCP<DataTransferKit::EntityLocalMap> local_map =
	Teuchos::rcp( new DataTransferKit::STKMeshEntityLocalMap(bulk_data) );

    // Create a DTK entity for the hex.
    DataTransferKit::Entity dtk_entity = 
	DataTransferKit::STKMeshEntity( hex_entity, bulk_data.ptr() );

    // Test the 
    TEST_EQUALITY( local_map->measure(dtk_entity), 8.0 );

    // Test the centroid.
    Teuchos::Array<double> centroid( space_dim, 0.0 );
    local_map->centroid( dtk_entity, centroid() );
    TEST_EQUALITY( centroid[0], 1.0 );
    TEST_EQUALITY( centroid[1], 1.0 );
    TEST_EQUALITY( centroid[2], 1.0 );

    // Make a good point and a bad point.
    Teuchos::Array<double> good_point( space_dim );
    good_point[0] = 0.5;
    good_point[1] = 1.5;
    good_point[2] = 1.0;
    Teuchos::Array<double> bad_point( space_dim );
    bad_point[0] = 0.75;
    bad_point[1] = -1.75;
    bad_point[2] = 0.35;

    // Test the reference frame safeguard.
    TEST_ASSERT(
	local_map->isSafeToMapToReferenceFrame(dtk_entity,good_point()) );
    TEST_ASSERT(
	!local_map->isSafeToMapToReferenceFrame(dtk_entity,bad_point()) );

    // Test the mapping to reference frame.
    Teuchos::Array<double> ref_good_point( space_dim );
    bool good_map = local_map->mapToReferenceFrame( 
	dtk_entity, good_point(), ref_good_point() );
    TEST_ASSERT( good_map );
    TEST_EQUALITY( ref_good_point[0], -0.5 );
    TEST_EQUALITY( ref_good_point[1], 0.5 );
    TEST_EQUALITY( ref_good_point[2], 0.0 );
			    
    Teuchos::Array<double> ref_bad_point( space_dim );
    bool bad_map = local_map->mapToReferenceFrame( 
	dtk_entity, bad_point(), ref_bad_point() );
    TEST_ASSERT( bad_map );

    // Test the point inclusion.
    TEST_ASSERT( local_map->checkPointInclusion(dtk_entity,ref_good_point()) );
    TEST_ASSERT( !local_map->checkPointInclusion(dtk_entity,ref_bad_point()) );

    // Test the map to physical frame.
    Teuchos::Array<double> phy_good_point( space_dim );
    local_map->mapToPhysicalFrame(dtk_entity,ref_good_point(),phy_good_point());
    TEST_EQUALITY( good_point[0], phy_good_point[0] );
    TEST_EQUALITY( good_point[1], phy_good_point[1] );
    TEST_EQUALITY( good_point[2], phy_good_point[2] );

    Teuchos::Array<double> phy_bad_point( space_dim );
    local_map->mapToPhysicalFrame(dtk_entity,ref_bad_point(),phy_bad_point());
    TEST_EQUALITY( bad_point[0], phy_bad_point[0] );
    TEST_EQUALITY( bad_point[1], phy_bad_point[1] );
    TEST_EQUALITY( bad_point[2], phy_bad_point[2] );
}

//---------------------------------------------------------------------------//
// end tstSTKMeshEntityLocalMap.cpp
//---------------------------------------------------------------------------//

