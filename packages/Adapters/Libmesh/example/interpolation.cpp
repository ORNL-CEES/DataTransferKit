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
 * \brief  Libmesh file-based interpolation example. Requires Libmesh to be
 * compiled with NetCDF reader.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "DTK_LibmeshManager.hpp"

#include "DTK_MapOperatorFactory.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <Tpetra_MultiVector.hpp>

#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/explicit_system.h>
#include <libmesh/fe.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/node.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/parallel.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/system.h>

//---------------------------------------------------------------------------//
// Data field function.
//---------------------------------------------------------------------------//
double dataFunction( double x, double y, double z )
{
    return std::abs( x ) + std::abs( y ) + std::abs( z ) + 1.0;
}

//---------------------------------------------------------------------------//
// Libmesh assemble functions.
//---------------------------------------------------------------------------//
void src_assemble_function( libMesh::EquationSystems &, const std::string & )
{ /* ... */}

void tgt_assemble_function( libMesh::EquationSystems &, const std::string & )
{ /* ... */}

//---------------------------------------------------------------------------//
// Example driver.
//---------------------------------------------------------------------------//
int main( int argc, char *argv[] )
{
    // INITIALIZATION
    // --------------

    // Setup communication.
    Teuchos::GlobalMPISession mpiSession( &argc, &argv );

    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();

    // Read in command line options.
    std::string xml_input_filename;
    Teuchos::CommandLineProcessor clp( false );
    clp.setOption( "xml-in-file", &xml_input_filename,
                   "The XML file to read into a parameter list" );
    clp.parse( argc, argv );

    // Build the parameter list from the xml input.
    Teuchos::RCP<Teuchos::ParameterList> plist =
        Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( xml_input_filename,
                                          Teuchos::inoutArg( *plist ) );

    // Read command-line options
    std::string src_mesh_input_file =
        plist->get<std::string>( "Source Mesh Input File" );
    std::string src_mesh_output_file =
        plist->get<std::string>( "Source Mesh Output File" );
    std::string tgt_mesh_input_file =
        plist->get<std::string>( "Target Mesh Input File" );
    std::string tgt_mesh_output_file =
        plist->get<std::string>( "Target Mesh Output File" );

    // Get the raw mpi communicator.
    Teuchos::RCP<const Teuchos::MpiComm<int>> mpi_comm =
        Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>( comm );
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> opaque_comm =
        mpi_comm->getRawMpiComm();
    MPI_Comm raw_comm = ( *opaque_comm )();

    // Initialize libmesh.
    libMesh::LibMeshInit libmesh_init( argc, argv, raw_comm );

    // SOURCE MESH READ
    // ----------------

    // Create a source mesh.
    Teuchos::RCP<libMesh::Mesh> src_mesh =
        Teuchos::rcp( new libMesh::Mesh( libmesh_init.comm() ) );
    src_mesh->read( src_mesh_input_file );

    // Make a libmesh system. We will put a first order linear basis on the
    // elements for all subdomains.
    std::string src_var_name = "src_var";
    libMesh::EquationSystems src_equation_systems( *src_mesh );
    libMesh::System &src_system =
        src_equation_systems.add_system<libMesh::ExplicitSystem>(
            "Source System" );
    int src_var_id = src_system.add_variable( src_var_name );
    src_system.attach_assemble_function( src_assemble_function );
    src_equation_systems.init();

    // Create some data and put it on the variable.
    libMesh::Mesh::const_element_iterator elements_begin =
        src_mesh->elements_begin();
    libMesh::Mesh::const_element_iterator elements_end =
        src_mesh->elements_end();
    int elem_n_nodes = 0;
    libMesh::Node *node;
    std::vector<libMesh::dof_id_type> elem_dof_ids;
    for ( auto element = elements_begin; element != elements_end; ++element )
    {
        src_system.get_dof_map().dof_indices( *element, elem_dof_ids,
                                              src_var_id );
        elem_n_nodes = ( *element )->n_nodes();

        for ( int n = 0; n < elem_n_nodes; ++n )
        {
            node = ( *element )->get_node( n );
            src_system.solution->set(
                elem_dof_ids[n], dataFunction( ( *node )( 0 ), ( *node )( 1 ),
                                               ( *node )( 2 ) ) );
        }
    }

    // TARGET MESH READ
    // ----------------

    // Create a target mesh.
    Teuchos::RCP<libMesh::Mesh> tgt_mesh =
        Teuchos::rcp( new libMesh::Mesh( libmesh_init.comm() ) );
    tgt_mesh->read( tgt_mesh_input_file );

    // Make a libmesh system. We will put a first order linear basis on the
    // elements for all subdomains.
    std::string tgt_var_name = "tgt_var";
    std::string tgt_err_name = "tgt_err";
    libMesh::EquationSystems tgt_equation_systems( *tgt_mesh );
    libMesh::System &tgt_system =
        tgt_equation_systems.add_system<libMesh::ExplicitSystem>(
            "Target System" );
    int tgt_var_id = tgt_system.add_variable( tgt_var_name );
    int tgt_err_id = tgt_system.add_variable( tgt_err_name );
    tgt_system.attach_assemble_function( tgt_assemble_function );
    tgt_equation_systems.init();

    // SOLUTION TRANSFER SETUP
    // -----------------------

    // Create a manager for the source set elements.
    DataTransferKit::LibmeshManager src_manager(
        src_mesh, Teuchos::rcpFromRef( src_system ) );

    // Create a manager for the target set nodes.
    DataTransferKit::LibmeshManager tgt_manager(
        tgt_mesh, Teuchos::rcpFromRef( tgt_system ) );

    // Create a solution vector for the source.
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        src_vector = src_manager.createFieldMultiVector( src_var_name );

    // Create a solution vector for the target.
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        tgt_vector = tgt_manager.createFieldMultiVector( tgt_var_name );

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

    // Create a map operator.
    Teuchos::ParameterList &dtk_list = plist->sublist( "DataTransferKit" );
    DataTransferKit::MapOperatorFactory op_factory;
    Teuchos::RCP<DataTransferKit::MapOperator> map_op = op_factory.create(
        src_vector->getMap(), tgt_vector->getMap(), dtk_list );

    // Setup the map operator.
    map_op->setup( src_manager.functionSpace(), tgt_manager.functionSpace() );

    // Apply the map operator.
    map_op->apply( *src_vector, *tgt_vector );

    // COMPUTE THE SOLUTION ERROR
    // --------------------------

    double gold_value = 0.0;
    double error_l2_norm = 0.0;
    double var_l2_norm = 0.0;
    double tgt_value = 0.0;

    std::unordered_map<libMesh::dof_id_type, double> tgt_val_map;
    std::unordered_map<libMesh::dof_id_type, double> error_val_map;
    elements_begin = tgt_mesh->elements_begin();
    elements_end = tgt_mesh->elements_end();
    elem_n_nodes = 0;
    std::vector<libMesh::dof_id_type> err_dof_ids;
    for ( auto element = elements_begin; element != elements_end; ++element )
    {
        tgt_system.get_dof_map().dof_indices( *element, elem_dof_ids,
                                              tgt_var_id );
        tgt_system.get_dof_map().dof_indices( *element, err_dof_ids,
                                              tgt_err_id );
        elem_n_nodes = ( *element )->n_nodes();

        for ( int n = 0; n < elem_n_nodes; ++n )
        {
            node = ( *element )->get_node( n );
            gold_value =
                dataFunction( ( *node )( 0 ), ( *node )( 1 ), ( *node )( 2 ) );
            tgt_value = tgt_system.solution->el( elem_dof_ids[n] );
            tgt_system.solution->set( err_dof_ids[n],
                                      ( tgt_value - gold_value ) / gold_value );
            tgt_val_map.emplace( err_dof_ids[n], tgt_value * tgt_value );
            error_val_map.emplace( err_dof_ids[n],
                                   ( tgt_value - gold_value ) *
                                       ( tgt_value - gold_value ) );
        }
    }

    for ( auto tval = tgt_val_map.begin(), eval = error_val_map.begin();
          tval != tgt_val_map.end(); ++tval, ++eval )
    {
        var_l2_norm += tval->second;
        error_l2_norm += eval->second;
    }

    error_l2_norm = std::sqrt( error_l2_norm );
    var_l2_norm = std::sqrt( var_l2_norm );
    std::cout << "|e|_2 / |f|_2: " << error_l2_norm / var_l2_norm << std::endl;

    // SOURCE MESH WRITE
    // -----------------

    libMesh::ExodusII_IO( *src_mesh )
        .write_equation_systems( src_mesh_output_file, src_equation_systems );

    // TARGET MESH WRITE
    // -----------------

    libMesh::ExodusII_IO( *tgt_mesh )
        .write_equation_systems( tgt_mesh_output_file, tgt_equation_systems );
}

//---------------------------------------------------------------------------//
// end interpolation.cpp
//---------------------------------------------------------------------------//
