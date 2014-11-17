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

    // Make a data field.
    stk::mesh::Field<double, stk::mesh::Cartesian3d>& data_field =
	meta_data.declare_field<
	stk::mesh::Field<double, stk::mesh::Cartesian3d> >(
	    stk::topology::NODE_RANK, "test field");
    meta_data.set_coordinate_field( &data_field );
    stk::mesh::put_field( data_field, part_1 );
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
    stk::mesh::Field<double,stk::mesh::Cartesian3d>* test_field =
	bulk_data->mesh_meta_data(
	    ).get_field<stk::mesh::Field<double,stk::mesh::Cartesian3d> >(
		stk::topology::NODE_RANK, "test field" );
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > field_vec =
	DataTransferKit::STKMeshDOFVector::createTpetraMultiVectorFromSTKField<double>(
	    *bulk_data, *test_field, 3 );

    // Test the vector.
    unsigned comm_size = comm->getSize();
    TEST_EQUALITY( 3, field_vec->getNumVectors() );
    TEST_EQUALITY( 8, field_vec->getLocalLength() );
    TEST_EQUALITY( 8*comm_size, field_vec->getGlobalLength() );
    
    // Put some data in the vector.
    double val_0 = 3.3;
    double val_1 = -9.3;
    double val_2 = 1.74;
    field_vec->getVectorNonConst( 0 )->putScalar( val_0 );
    field_vec->getVectorNonConst( 1 )->putScalar( val_1 );
    field_vec->getVectorNonConst( 2 )->putScalar( val_2 );

    // Push the data back to STK.
    DataTransferKit::STKMeshDOFVector::pushTpetraMultiVectorToSTKField(
    	field_vec, *bulk_data, *test_field );

    // Test the STK field.
    for ( stk::mesh::Entity node : nodes )
    {
    	double* data = stk::mesh::field_data( *test_field, node );
    	TEST_EQUALITY( data[0], val_0 );
    	TEST_EQUALITY( data[1], val_1 );
    	TEST_EQUALITY( data[2], val_2 );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( STKMeshDOFVector, view_test )
{
    // Initialize parallel communication.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();

    // Vector parameters.
    int num_vec = 3;
    int vec_length = 10;

    // Create data.
    Teuchos::Array<std::size_t> ids( vec_length );
    Teuchos::ArrayRCP<double> in_data( num_vec*vec_length );
    Teuchos::ArrayRCP<double> out_data( num_vec*vec_length );
    for ( int i = 0; i < vec_length; ++i )
    {
	ids[i] = i+1;
	in_data[i] = 2.0*i + 1.0;
	out_data[i] = 0.0;
    }

    // Create an input vector.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > in_vec =
	DataTransferKit::STKMeshDOFVector::createTpetraMultiVectorFromView( 
	    comm, ids(), in_data, vec_length, num_vec );
        
    // Create an output vector.
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > out_vec =
	DataTransferKit::STKMeshDOFVector::createTpetraMultiVectorFromView( 
	    comm, ids(), out_data, vec_length, num_vec );

    // Add the vectors together.
    out_vec->update( 1.0, *in_vec, 0.0 );

    // Check the results.
    for ( int i = 0; i < vec_length; ++i )
    {
	TEST_EQUALITY( in_data[i], out_data[i] );
    }
}

//---------------------------------------------------------------------------//
// end of tstSTKMeshDOFVector.cpp
//---------------------------------------------------------------------------//
