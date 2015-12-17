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
 * \file   interpolation.cpp
 * \author Stuart Slattery
 * \brief  STK file-based interpolation example.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <cstdlib>

#include "DTK_STKMeshHelpers.hpp"
#include "DTK_STKMeshManager.hpp"
#include "DTK_MapOperatorFactory.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <Tpetra_MultiVector.hpp>

#include <Intrepid_FieldContainer.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/fixtures/HexFixture.hpp>
#include <stk_mesh/fixtures/CoordinateMapping.hpp>
#include <stk_topology/topology.hpp>

//---------------------------------------------------------------------------//
// Data field function.
//---------------------------------------------------------------------------//
double dataFunction( double x, double y, double z )
{
    return std::abs(x) + std::abs(y) + std::abs(z) + 1.0;
}

//---------------------------------------------------------------------------//
// Example driver.
//---------------------------------------------------------------------------//
int main(int argc, char* argv[])
{
    // INITIALIZATION
    // --------------

    // Setup communication.
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();

    // Get the raw mpi communicator (basic typedef in STK).
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm = 
	Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( comm );
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm = 
	mpi_comm->getRawMpiComm();
    stk::ParallelMachine parallel_machine = (*opaque_comm)();

    // Read in command line options.
    std::string xml_input_filename;
    Teuchos::CommandLineProcessor clp(false);
    clp.setOption( "xml-in-file",
		   &xml_input_filename,
		   "The XML file to read into a parameter list" );
    clp.parse(argc,argv);

    // Build the parameter list from the xml input.
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile(
	xml_input_filename, Teuchos::inoutArg(*plist) );

    // Create a bounding box for the test.
    int x_max = 10;
    int y_max = 10;
    int z_max = 10;
    int src_mesh_size = 10;
    int tgt_mesh_size = 11;
    stk::mesh::fixtures::FixedCartesianCoordinateMapping 
	src_coord_mapping( src_mesh_size, src_mesh_size, src_mesh_size, x_max, y_max, z_max );
    stk::mesh::fixtures::FixedCartesianCoordinateMapping 
	tgt_coord_mapping( tgt_mesh_size, tgt_mesh_size, tgt_mesh_size, x_max, y_max, z_max );
    
    // SOURCE MESH READ
    // ----------------

    // Create a box fixture.
    stk::mesh::fixtures::HexFixture src_fixture( parallel_machine,
						 src_mesh_size, src_mesh_size, src_mesh_size );

    // Add a nodal field to the source part.
    stk::mesh::Field<double>& source_field = 
	src_fixture.m_meta.declare_field<stk::mesh::Field<double> >( 
	    stk::topology::NODE_RANK, "u_src" );
    stk::mesh::Part& src_part = src_fixture.m_meta.universal_part();
    stk::mesh::put_field( source_field, src_part );

    // Create the mesh on the fixture.
    Teuchos::RCP<stk::mesh::BulkData> src_bulk_data = 
	Teuchos::rcpFromRef( src_fixture.m_bulk_data );
    src_bulk_data->modification_begin();
    src_fixture.generate_mesh( src_coord_mapping );
    src_bulk_data->modification_end();

    // Put some data in the source field. We will use the distance of the node
    // from the origin as the data.
    stk::mesh::Selector src_stk_selector( src_part );
    stk::mesh::BucketVector src_part_buckets = 
	src_stk_selector.get_buckets( stk::topology::NODE_RANK );
    std::vector<stk::mesh::Entity> src_part_nodes;
    stk::mesh::get_selected_entities( 
	src_stk_selector, src_part_buckets, src_part_nodes );
    Intrepid::FieldContainer<double> src_node_coords =
	DataTransferKit::STKMeshHelpers::getEntityNodeCoordinates(
	    Teuchos::Array<stk::mesh::Entity>(src_part_nodes), *src_bulk_data );
    double* src_field_data;
    int num_src_part_nodes = src_part_nodes.size();
    for ( int n = 0; n < num_src_part_nodes; ++n )
    {
	src_field_data = stk::mesh::field_data( source_field, src_part_nodes[n] );
	src_field_data[0] = dataFunction( src_node_coords(n,0,0),
					  src_node_coords(n,0,1),
					  src_node_coords(n,0,2) );
    }


    // TARGET MESH READ
    // ----------------

    // Create a box fixture.
    stk::mesh::fixtures::HexFixture tgt_fixture( parallel_machine,
						 tgt_mesh_size, tgt_mesh_size, tgt_mesh_size );

    // Get the target part.
    stk::mesh::Part& tgt_part = tgt_fixture.m_meta.universal_part();
    stk::mesh::Selector tgt_stk_selector( tgt_part );

    // Add a nodal field to the target part.
    stk::mesh::Field<double>& target_field = 
    	tgt_fixture.m_meta.declare_field<stk::mesh::Field<double> >( 
    	    stk::topology::NODE_RANK, "u_tgt" );
    stk::mesh::put_field( target_field, tgt_part );

    // Add an error nodal field to the target part.
    stk::mesh::Field<double>& target_error_field = 
    	tgt_fixture.m_meta.declare_field<stk::mesh::Field<double> >( 
    	    stk::topology::NODE_RANK, "u_err" );
    stk::mesh::put_field( target_error_field, tgt_part );

    // Create the mesh on the fixture.
    Teuchos::RCP<stk::mesh::BulkData> tgt_bulk_data = 
	Teuchos::rcpFromRef( tgt_fixture.m_bulk_data );
    tgt_bulk_data->modification_begin();
    tgt_fixture.generate_mesh( tgt_coord_mapping );
    tgt_bulk_data->modification_end();

    
    // SOLUTION TRANSFER SETUP
    // -----------------------
    
    // Create a manager for the source part elements.
    DataTransferKit::STKMeshManager src_manager( src_bulk_data, src_stk_selector );

    // Create a manager for the target part nodes.
    DataTransferKit::STKMeshManager tgt_manager( tgt_bulk_data, tgt_stk_selector );

    // Create a solution vector for the source.
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > src_vector =
	src_manager.createFieldMultiVector<stk::mesh::Field<double> >(
	    Teuchos::ptr(&source_field), 1 );
    
    // Create a solution vector for the target.
    Teuchos::RCP<Tpetra::MultiVector<double,int,DataTransferKit::SupportId> > tgt_vector =
	tgt_manager.createFieldMultiVector<stk::mesh::Field<double> >(
	    Teuchos::ptr(&target_field), 1 );

    // Print out source mesh info.
    Teuchos::RCP<Teuchos::Describable> src_describe =
	src_manager.functionSpace()->entitySet();
    std::cout << "Source Mesh" << std::endl;
    src_describe->describe( std::cout );
    std::cout << std::endl;

    // Print out target mesh info.
    Teuchos::RCP<Teuchos::Describable> tgt_describe =
	tgt_manager.functionSpace()->entitySet();
    std::cout << "Target Mesh" << std::endl;
    tgt_describe->describe( std::cout );
    std::cout << std::endl;

    
    // SOLUTION TRANSFER
    // -----------------

    // Create a map operator. The operator settings are in the
    // "DataTransferKit" parameter list.
    Teuchos::ParameterList& dtk_list = plist->sublist("DataTransferKit");    
    DataTransferKit::MapOperatorFactory op_factory;
    Teuchos::RCP<DataTransferKit::MapOperator> map_op =
	op_factory.create( src_vector->getMap(),
			   tgt_vector->getMap(),
			   dtk_list );

    // Setup the map operator. This creates the underlying linear operators.
    map_op->setup( src_manager.functionSpace(), tgt_manager.functionSpace() );

    // Apply the map operator. This interpolates the data from one STK field
    // to the other.
    map_op->apply( *src_vector, *tgt_vector );


    // COMPUTE THE SOLUTION ERROR
    // --------------------------

    stk::mesh::BucketVector tgt_part_buckets = 
	tgt_stk_selector.get_buckets( stk::topology::NODE_RANK );
    std::vector<stk::mesh::Entity> tgt_part_nodes;
    stk::mesh::get_selected_entities( 
	tgt_stk_selector, tgt_part_buckets, tgt_part_nodes );
    Intrepid::FieldContainer<double> tgt_node_coords =
	DataTransferKit::STKMeshHelpers::getEntityNodeCoordinates(
	    Teuchos::Array<stk::mesh::Entity>(tgt_part_nodes), *tgt_bulk_data );
    double* tgt_field_data;
    double* err_field_data;
    int num_tgt_part_nodes = tgt_part_nodes.size();
    double error_l2_norm = 0.0;
    double field_l2_norm = 0.0;
    for ( int n = 0; n < num_tgt_part_nodes; ++n )
    {
	double gold_value = dataFunction( tgt_node_coords(n,0,0),
					  tgt_node_coords(n,0,1),
					  tgt_node_coords(n,0,2) );
	tgt_field_data = stk::mesh::field_data( target_field, tgt_part_nodes[n] );
	err_field_data = stk::mesh::field_data( target_error_field, tgt_part_nodes[n] );
	err_field_data[0] = tgt_field_data[0] - gold_value;
	error_l2_norm += err_field_data[0] * err_field_data[0];
	field_l2_norm += tgt_field_data[0] * tgt_field_data[0];
	err_field_data[0] /= gold_value;
    }
    error_l2_norm = std::sqrt( error_l2_norm );
    field_l2_norm = std::sqrt( field_l2_norm );
    std::cout << "|e|_2 / |f|_2: " << error_l2_norm / field_l2_norm << std::endl;
}

//---------------------------------------------------------------------------//
// end tstSTK_Mesh.cpp
//---------------------------------------------------------------------------//
