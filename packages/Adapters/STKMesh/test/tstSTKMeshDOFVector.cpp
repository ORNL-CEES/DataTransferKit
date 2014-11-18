//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tstSTKMeshDOFVector.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  Entity-centered DOF vector test.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "DTK_STKMeshDOFVector.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include <Teuchos_DefaultMpiComm.hpp>

#include <Tpetra_MultiVector.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_topology/topology.hpp>

//---------------------------------------------------------------------------//
// Hex-8 test.
TEUCHOS_UNIT_TEST( STKMeshEntitySet, pull_push_test )
{
    // Extract the raw mpi communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
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

    std::string p2_name = "part_2";
    stk::mesh::Part& part_2 = meta_data.declare_part( p2_name );

    // Make a data field.
    stk::mesh::Field<double, stk::mesh::Cartesian3d>& data_field_1 =
	meta_data.declare_field<
	stk::mesh::Field<double, stk::mesh::Cartesian3d> >(
	    stk::topology::NODE_RANK, "test field 1");
    meta_data.set_coordinate_field( &data_field_1 );
    stk::mesh::put_field( data_field_1, part_1 );

    // Make an empty data field.
    stk::mesh::Field<double, stk::mesh::Cartesian3d>& data_field_2 =
	meta_data.declare_field<
	stk::mesh::Field<double, stk::mesh::Cartesian3d> >(
	    stk::topology::NODE_RANK, "test field 2");
    meta_data.set_coordinate_field( &data_field_2 );
    stk::mesh::put_field( data_field_2, part_2 );

    // Commit the meta data.
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

    // Create a vector from the nodal field.
    stk::mesh::Field<double,stk::mesh::Cartesian3d>* test_field_1 =
	bulk_data->mesh_meta_data(
	    ).get_field<stk::mesh::Field<double,stk::mesh::Cartesian3d> >(
		stk::topology::NODE_RANK, "test field 1" );
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > field_vec_1 =
	DataTransferKit::STKMeshDOFVector::createTpetraMultiVectorFromSTKField<double>(
	    *bulk_data, *test_field_1, 3 );

    // Test the vector.
    unsigned comm_size = comm->getSize();
    TEST_EQUALITY( 3, field_vec_1->getNumVectors() );
    TEST_EQUALITY( 8, field_vec_1->getLocalLength() );
    TEST_EQUALITY( 8*comm_size, field_vec_1->getGlobalLength() );
    
    // Put some data in the vector.
    double val_0 = 3.3;
    double val_1 = -9.3;
    double val_2 = 1.74;
    field_vec_1->getVectorNonConst( 0 )->putScalar( val_0 );
    field_vec_1->getVectorNonConst( 1 )->putScalar( val_1 );
    field_vec_1->getVectorNonConst( 2 )->putScalar( val_2 );

    // Push the data back to STK.
    DataTransferKit::STKMeshDOFVector::pushTpetraMultiVectorToSTKField(
    	*field_vec_1, *bulk_data, *test_field_1 );

    // Test the STK field.
    for ( stk::mesh::Entity node : nodes )
    {
    	double* data = stk::mesh::field_data( *test_field_1, node );
    	TEST_EQUALITY( data[0], val_0 );
    	TEST_EQUALITY( data[1], val_1 );
    	TEST_EQUALITY( data[2], val_2 );
    }

    // Now make an empty field vector.
    stk::mesh::Field<double,stk::mesh::Cartesian3d>* test_field_2 =
	bulk_data->mesh_meta_data(
	    ).get_field<stk::mesh::Field<double,stk::mesh::Cartesian3d> >(
		stk::topology::NODE_RANK, "test field 2" );
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > field_vec_2 =
	DataTransferKit::STKMeshDOFVector::createTpetraMultiVectorFromSTKField<double>(
	    *bulk_data, *test_field_2, 3 );

    // Test the vector to make sure it is empty.
    TEST_EQUALITY( 3, field_vec_2->getNumVectors() );
    TEST_EQUALITY( 0, field_vec_2->getLocalLength() );
    TEST_EQUALITY( 0, field_vec_2->getGlobalLength() );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( STKMeshDOFVector, view_test )
{
    // Extract the raw mpi communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
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

    // Commit the meta data.
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

    // Get the nodes.
    std::vector<stk::mesh::Entity> mesh_nodes;
    stk::mesh::get_entities( *bulk_data, stk::topology::NODE_RANK, mesh_nodes );
    TEST_EQUALITY( num_nodes, mesh_nodes.size() );

    // Vector parameters.
    int field_dim = 3;

    // Create data.
    Teuchos::Array<std::size_t> ids( num_nodes );
    Teuchos::ArrayRCP<double> in_data( field_dim*num_nodes );
    Teuchos::ArrayRCP<double> out_data( field_dim*num_nodes );
    for ( unsigned i = 0; i < num_nodes; ++i )
    {
	ids[i] = i+1;
	in_data[i] = 2.0*i + 1.0;
	out_data[i] = 0.0;
    }

    // Create a vector from the entities and view.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > vec_1 =
	DataTransferKit::STKMeshDOFVector::createTpetraMultiVectorFromEntitiesAndView( 
	    *bulk_data, mesh_nodes, field_dim, in_data );
    
    // Test the vector.
    unsigned comm_size = comm->getSize();
    TEST_EQUALITY( 3, vec_1->getNumVectors() );
    TEST_EQUALITY( num_nodes, vec_1->getLocalLength() );
    TEST_EQUALITY( num_nodes*comm_size, vec_1->getGlobalLength() );
        
    // Create a second vector from the entities and view.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > vec_2 =
	DataTransferKit::STKMeshDOFVector::createTpetraMultiVectorFromEntitiesAndView( 
	    *bulk_data, mesh_nodes, field_dim, out_data );

    // Test the vector.
    TEST_EQUALITY( 3, vec_2->getNumVectors() );
    TEST_EQUALITY( num_nodes, vec_2->getLocalLength() );
    TEST_EQUALITY( num_nodes*comm_size, vec_2->getGlobalLength() );

    // Add vector 1 and vector 2 together.
    vec_2->update( 1.0, *vec_1, 0.0 );

    // Check the results.
    for ( unsigned i = 0; i < num_nodes; ++i )
    {
	TEST_EQUALITY( in_data[i], out_data[i] );
    }

    // Create another Tpetra vector from a part vector.
    stk::mesh::PartVector parts( 1, &part_1 );
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > vec_3 =
	DataTransferKit::STKMeshDOFVector::createTpetraMultiVectorFromPartVectorAndView( 
	    *bulk_data, parts, stk::topology::NODE_RANK, field_dim, out_data );

    // Test the vector.
    TEST_EQUALITY( 3, vec_3->getNumVectors() );
    TEST_EQUALITY( num_nodes, vec_3->getLocalLength() );
    TEST_EQUALITY( num_nodes*comm_size, vec_3->getGlobalLength() );
    vec_3->putScalar( 9.9 );
    for ( unsigned i = 0; i < num_nodes; ++i )
    {
	TEST_EQUALITY( 9.9, out_data[i] );
    }

    // Create another Tpetra vector from a selector.
    stk::mesh::Selector selector( part_1 );
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > vec_4 =
	DataTransferKit::STKMeshDOFVector::createTpetraMultiVectorFromSelectorAndView( 
	    *bulk_data, selector, stk::topology::NODE_RANK, field_dim, in_data );

    // Test the vector.
    TEST_EQUALITY( 3, vec_4->getNumVectors() );
    TEST_EQUALITY( num_nodes, vec_4->getLocalLength() );
    TEST_EQUALITY( num_nodes*comm_size, vec_4->getGlobalLength() );

    // Add vector 3 and vector 4 together.
    vec_4->update( 1.0, *vec_3, 0.0 );

    // Check the results.
    for ( unsigned i = 0; i < num_nodes; ++i )
    {
	TEST_EQUALITY( in_data[i], out_data[i] );
    }
}

//---------------------------------------------------------------------------//
// end of tstSTKMeshDOFVector.cpp
//---------------------------------------------------------------------------//
