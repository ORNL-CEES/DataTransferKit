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
 * \file   tstMoabTagField.cpp
 * \author Stuart Slattery
 * \brief  Tag field test
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_MoabTagField.hpp>
#include <DTK_FieldMultiVector.hpp>
#include <DTK_MoabEntitySet.hpp>
#include <DTK_MoabMeshSetIndexer.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

#include <Tpetra_MultiVector.hpp>

#include <moab/Interface.hpp>
#include <moab/ParallelComm.hpp>
#include <moab/Core.hpp>

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( MoabTagField, push_pull_test )
{
    // Extract the raw mpi communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm = 
	Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm<int> >( comm );
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm = 
	mpi_comm->getRawMpiComm();
    MPI_Comm raw_comm = (*opaque_comm)();

    // Create the mesh.
    Teuchos::RCP<moab::Interface> moab_mesh = Teuchos::rcp( new moab::Core() );
    Teuchos::RCP<moab::ParallelComm> parallel_mesh =
    	Teuchos::rcp( new moab::ParallelComm(moab_mesh.getRawPtr(),raw_comm) );

    // Create the nodes.
    moab::ErrorCode error = moab::MB_SUCCESS;
    unsigned num_nodes = 8;
    Teuchos::Array<moab::EntityHandle> nodes(num_nodes);
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

    // Make 2 entity sets.
    moab::EntityHandle entity_set_1;
    error = moab_mesh->create_meshset( 0, entity_set_1 );
    TEST_EQUALITY( error, moab::MB_SUCCESS );
    moab::EntityHandle entity_set_2;
    error = moab_mesh->create_meshset( 0, entity_set_2 );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    // Put the nodes in the first entity set.
    error = moab_mesh->add_entities( entity_set_1, nodes.getRawPtr(), num_nodes );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    // Create a tag for entity set 1.
    int tag_size = 3;
    moab::Tag tag_1;
    bool created = false;
    double default_val = 0.0;
    Teuchos::Array<double> default_tag( tag_size, default_val );
    error = moab_mesh->tag_get_handle( 
    	"Tag_1",
    	tag_size,
    	moab::MB_TYPE_DOUBLE,
    	tag_1,
    	moab::MB_TAG_DENSE|moab::MB_TAG_CREAT,
	static_cast<void*>(default_tag.getRawPtr()),
	&created );
    TEST_EQUALITY( error, moab::MB_SUCCESS );
    TEST_ASSERT( created );

    double test_val = 2.2;
    Teuchos::Array<double> test_tag( num_nodes*tag_size, test_val );
    error = moab_mesh->tag_set_data( 
    	tag_1,
    	nodes.getRawPtr(),
    	num_nodes,
    	static_cast<void*>(test_tag.getRawPtr()) );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    // Create an entity set.
    Teuchos::RCP<DataTransferKit::MoabMeshSetIndexer> set_indexer = Teuchos::rcp(
	new DataTransferKit::MoabMeshSetIndexer(parallel_mesh) );
    Teuchos::RCP<DataTransferKit::EntitySet> dtk_entity_set = Teuchos::rcp(
	new DataTransferKit::MoabEntitySet(parallel_mesh,set_indexer) );

    // Create a vector from entity set 1.
    Teuchos::RCP<DataTransferKit::Field> field_1 = Teuchos::rcp(
	new DataTransferKit::MoabTagField<double>(
	    parallel_mesh, set_indexer, entity_set_1, tag_1) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > tag_vec_1 =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  field_1, dtk_entity_set) );

    // Test the vector size.
    unsigned comm_size = comm->getSize();
    TEST_EQUALITY( tag_size, Teuchos::as<int>(tag_vec_1->getNumVectors()) );
    TEST_EQUALITY( num_nodes, tag_vec_1->getLocalLength() );
    TEST_EQUALITY( num_nodes*comm_size, tag_vec_1->getGlobalLength() );

    // Test the default tag data.
    Teuchos::rcp_dynamic_cast<DataTransferKit::FieldMultiVector>(
	tag_vec_1)->pullDataFromApplication();
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double> > tag_vec_1_view =
    		      tag_vec_1->get2dView();
    for ( unsigned n = 0; n < num_nodes; ++n )
    {
    	TEST_EQUALITY( tag_vec_1_view[0][n], test_val );
    	TEST_EQUALITY( tag_vec_1_view[1][n], test_val );
    	TEST_EQUALITY( tag_vec_1_view[2][n], test_val );
    }
    
    // Put some data in the vector.
    double val_0 = 3.3;
    double val_1 = -9.3;
    double val_2 = 1.74;
    tag_vec_1->getVectorNonConst( 0 )->putScalar( val_0 );
    tag_vec_1->getVectorNonConst( 1 )->putScalar( val_1 );
    tag_vec_1->getVectorNonConst( 2 )->putScalar( val_2 );

    // Push the data back to Moab
    Teuchos::rcp_dynamic_cast<DataTransferKit::FieldMultiVector>(
	tag_vec_1)->pushDataToApplication();

    // Test the Moab tag.
    Teuchos::Array<const void*> node_data( num_nodes );
    error = moab_mesh->tag_get_by_ptr( tag_1,
    				       nodes.getRawPtr(),
    				       num_nodes,
    				       node_data.getRawPtr() );
    for ( unsigned n = 0; n < num_nodes; ++n )
    {
    	const double* data = static_cast<const double*>(node_data[n]);
    	TEST_EQUALITY( data[0], val_0 );
    	TEST_EQUALITY( data[1], val_1 );
    	TEST_EQUALITY( data[2], val_2 );
    }

    // Create a tag for entity set 2.
    moab::Tag tag_2;
    error = moab_mesh->tag_get_handle( 
    	"Tag_2",
    	tag_size,
    	moab::MB_TYPE_DOUBLE,
    	tag_2,
    	moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    // Make an empty vector over set 2.
    Teuchos::RCP<DataTransferKit::Field> field_2 = Teuchos::rcp(
	new DataTransferKit::MoabTagField<double>(
	    parallel_mesh, set_indexer, entity_set_2, tag_2) );
    Teuchos::RCP<Tpetra::MultiVector<double,int,std::size_t> > tag_vec_2 =
	Teuchos::rcp( new DataTransferKit::FieldMultiVector(
			  field_2, dtk_entity_set) );

    // Test the vector to make sure it is empty.
    TEST_EQUALITY( tag_size, Teuchos::as<int>(tag_vec_2->getNumVectors()) );
    TEST_EQUALITY( 0, tag_vec_2->getLocalLength() );
    TEST_EQUALITY( 0, tag_vec_2->getGlobalLength() );
}

//---------------------------------------------------------------------------//
// end of tstMoabTagField.cpp
//---------------------------------------------------------------------------//
