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
//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tstSTKMeshField.cpp
 * \author Stuart Slattery
 * \brief  Field test.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "DTK_FieldMultiVector.hpp"
#include "DTK_STKMeshManager.hpp"

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
TEUCHOS_UNIT_TEST( STKMeshField, pull_push_test )
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

    // Create a nodal field.
    stk::mesh::Field<double,stk::mesh::Cartesian3d>* test_field_1 =
	bulk_data->mesh_meta_data(
	    ).get_field<stk::mesh::Field<double,stk::mesh::Cartesian3d> >(
		stk::topology::NODE_RANK, "test field 1" );
    for ( stk::mesh::Entity node : nodes )
    {
    	double* data = stk::mesh::field_data( *test_field_1, node );
	data[0] = 1.0;
	data[1] = 2.0;
	data[2] = 3.0;
    }

    // Create a manager.
    DataTransferKit::STKMeshManager manager( bulk_data );
    
    // Create a vector from the nodal field.
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > field_vec_1 =
	manager.createFieldMultiVector<stk::mesh::Field<double,stk::mesh::Cartesian3d> >(
	    Teuchos::ptr(test_field_1), 3 );

    // Test the vector allocation.
    unsigned comm_size = comm->getSize();
    TEST_EQUALITY( 3, field_vec_1->getNumVectors() );
    TEST_EQUALITY( 8, field_vec_1->getLocalLength() );
    TEST_EQUALITY( 8*comm_size, field_vec_1->getGlobalLength() );

    // Test the vector data.
    Teuchos::rcp_dynamic_cast<DataTransferKit::FieldMultiVector>(
	field_vec_1)->pullDataFromApplication();
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double> > field_vec_1_view =
		      field_vec_1->get2dView();
    for ( unsigned n = 0; n < num_nodes; ++n )
    {
	TEST_EQUALITY( field_vec_1_view[0][n], 1.0 );
	TEST_EQUALITY( field_vec_1_view[1][n], 2.0 );
	TEST_EQUALITY( field_vec_1_view[2][n], 3.0 );
    }
    
    // Put some data in the vector.
    double val_0 = 3.3;
    double val_1 = -9.3;
    double val_2 = 1.74;
    field_vec_1->getVectorNonConst( 0 )->putScalar( val_0 );
    field_vec_1->getVectorNonConst( 1 )->putScalar( val_1 );
    field_vec_1->getVectorNonConst( 2 )->putScalar( val_2 );

    // Push the data back to STK.
    Teuchos::rcp_dynamic_cast<DataTransferKit::FieldMultiVector>(
	field_vec_1)->pushDataToApplication();

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

    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > field_vec_2 =
	manager.createFieldMultiVector<stk::mesh::Field<double,stk::mesh::Cartesian3d> >(
	    Teuchos::ptr(test_field_2), 3 );

    // Test the vector to make sure it is empty.
    TEST_EQUALITY( 3, field_vec_2->getNumVectors() );
    TEST_EQUALITY( 0, field_vec_2->getLocalLength() );
    TEST_EQUALITY( 0, field_vec_2->getGlobalLength() );
}

//---------------------------------------------------------------------------//
// end of tstSTKMeshField.cpp
//---------------------------------------------------------------------------//
