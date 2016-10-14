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
 * \brief  Moab file-based interpolation example. Requires Moab to be compiled
 * with NetCDF reader.
 */
//---------------------------------------------------------------------------//

#include <moab/Core.hpp>
#include <moab/Interface.hpp>
#include <moab/ParallelComm.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <vector>

#include "DTK_MapOperatorFactory.hpp"
#include "DTK_MoabManager.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultComm.hpp>
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

//---------------------------------------------------------------------------//
// MOAB error check.
//---------------------------------------------------------------------------//
void checkMoabErrorCode( const moab::ErrorCode &ec )
{
    if ( moab::MB_SUCCESS != ec )
    {
        std::cout << "MOAB ERROR CODE " << ec << std::endl;
    }
}

//---------------------------------------------------------------------------//
// Data field function.
//---------------------------------------------------------------------------//
double dataFunction( double x, double y, double z )
{
    return std::abs( x ) + std::abs( y ) + std::abs( z ) + 1.0;
}

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
    std::string source_mesh_input_file =
        plist->get<std::string>( "Source Mesh Input File" );
    std::string source_mesh_output_file =
        plist->get<std::string>( "Source Mesh Output File" );
    std::string target_mesh_input_file =
        plist->get<std::string>( "Target Mesh Input File" );
    std::string target_mesh_output_file =
        plist->get<std::string>( "Target Mesh Output File" );

    // SOURCE MESH READ
    // ----------------

    // Create a source database.
    Teuchos::RCP<moab::Interface> source_iface =
        Teuchos::rcp( new moab::Core() );

    // Load the mesh.
    moab::ErrorCode error;
    error = source_iface->load_file( source_mesh_input_file.c_str(), 0, 0 );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );

    // Get the parallel moab instance.
    Teuchos::RCP<moab::ParallelComm> source_mesh = Teuchos::rcp(
        moab::ParallelComm::get_pcomm( source_iface.getRawPtr(), 0 ), false );

    // Get the entity set for the source part. Just use the root set for now.
    moab::EntityHandle source_set = source_iface->get_root_set();

    // Get the nodes in the source set.
    std::vector<moab::EntityHandle> source_nodes;
    error = source_iface->get_entities_by_type( source_set, moab::MBVERTEX,
                                                source_nodes );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );
    int num_source_nodes = source_nodes.size();

    // Create a mesh set for just the source nodes. We need this so we can
    // build a vector over the nodes.
    moab::EntityHandle source_node_set;
    error = source_iface->create_meshset( 0, source_node_set );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );
    error = source_iface->add_entities( source_node_set, source_nodes.data(),
                                        num_source_nodes );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );

    // Create a source data tag.
    moab::Tag source_data_tag;
    std::string source_data_tag_name( "u_src" );
    double default_val = 0.0;
    bool created = false;
    error = source_iface->tag_get_handle(
        source_data_tag_name.c_str(), 1, moab::MB_TYPE_DOUBLE, source_data_tag,
        moab::MB_TAG_DENSE | moab::MB_TAG_CREAT,
        static_cast<void *>( &default_val ), &created );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );
    assert( created );

    // Create some data and put it on the tag.
    Teuchos::Array<double> source_tag_data( num_source_nodes );
    Teuchos::Array<double> source_coords( 3 );
    for ( int n = 0; n < num_source_nodes; ++n )
    {
        error = source_iface->get_coords( &source_nodes[n], 1,
                                          source_coords.getRawPtr() );
        checkMoabErrorCode( error );
        assert( moab::MB_SUCCESS == error );
        source_tag_data[n] = dataFunction( source_coords[0], source_coords[1],
                                           source_coords[2] );
    }
    error = source_iface->tag_set_data(
        source_data_tag, source_nodes.data(), num_source_nodes,
        static_cast<void *>( source_tag_data.getRawPtr() ) );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );

    // TARGET MESH READ
    // ----------------

    // Create a target database.
    Teuchos::RCP<moab::Interface> target_iface =
        Teuchos::rcp( new moab::Core() );

    // Load the mesh.
    error = source_iface->load_file( target_mesh_input_file.c_str(), 0, 0 );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );

    // Get the parallel moab instance.
    Teuchos::RCP<moab::ParallelComm> target_mesh = Teuchos::rcp(
        moab::ParallelComm::get_pcomm( target_iface.getRawPtr(), 0 ), false );

    // Get the entity set for the target part. Just use the root set for now.
    moab::EntityHandle target_set = target_iface->get_root_set();

    // Get the nodes in the target set.
    std::vector<moab::EntityHandle> target_nodes;
    error = target_iface->get_entities_by_type( target_set, moab::MBVERTEX,
                                                target_nodes );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );
    int num_target_nodes = target_nodes.size();

    // Create a mesh set for just the target nodes. We need this so we can
    // build a vector over the nodes.
    moab::EntityHandle target_node_set;
    error = target_iface->create_meshset( 0, target_node_set );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );
    error = target_iface->add_entities( target_node_set, target_nodes.data(),
                                        num_target_nodes );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );

    // Create a target data tag.
    moab::Tag target_data_tag;
    std::string target_data_tag_name( "u_tgt" );
    error = target_iface->tag_get_handle(
        target_data_tag_name.c_str(), 1, moab::MB_TYPE_DOUBLE, target_data_tag,
        moab::MB_TAG_DENSE | moab::MB_TAG_CREAT,
        static_cast<void *>( &default_val ), &created );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );
    assert( created );

    // Create some data and put it on the tag.
    Teuchos::Array<double> target_tag_data( num_target_nodes, 0.0 );
    error = target_iface->tag_set_data(
        target_data_tag, target_nodes.data(), num_target_nodes,
        static_cast<void *>( target_tag_data.getRawPtr() ) );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );

    // Create a target error tag.
    moab::Tag target_error_tag;
    std::string target_error_tag_name( "u_err" );
    error = target_iface->tag_get_handle(
        target_error_tag_name.c_str(), 1, moab::MB_TYPE_DOUBLE,
        target_error_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT,
        static_cast<void *>( &default_val ), &created );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );
    assert( created );

    // SOLUTION TRANSFER SETUP
    // -----------------------

    // Create a manager for the source set elements.
    DataTransferKit::MoabManager src_manager( source_mesh, source_set );

    // Create a manager for the target set nodes.
    DataTransferKit::MoabManager tgt_manager( target_mesh, target_set );

    // Create a solution vector for the source.
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        src_vector = src_manager.createFieldMultiVector( source_node_set,
                                                         source_data_tag );

    // Create a solution vector for the target.
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        tgt_vector = tgt_manager.createFieldMultiVector( target_node_set,
                                                         target_data_tag );

    // Print out mesh info.
    Teuchos::RCP<Teuchos::FancyOStream> fancy_out =
        Teuchos::VerboseObjectBase::getDefaultOStream();
    fancy_out->setShowProcRank( true );
    src_manager.functionSpace()->entitySet()->describe( *fancy_out );
    tgt_manager.functionSpace()->entitySet()->describe( *fancy_out );

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
    double tag_l2_norm = 0.0;
    Teuchos::Array<double> error_tag_data( num_target_nodes );
    Teuchos::Array<double> target_coords( 3 );

    error = target_iface->tag_get_data(
        target_data_tag, target_nodes.data(), num_target_nodes,
        static_cast<void *>( target_tag_data.getRawPtr() ) );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );

    for ( int n = 0; n < num_target_nodes; ++n )
    {
        error = target_iface->get_coords( &target_nodes[n], 1,
                                          target_coords.getRawPtr() );
        checkMoabErrorCode( error );
        assert( moab::MB_SUCCESS == error );
        gold_value = dataFunction( target_coords[0], target_coords[1],
                                   target_coords[2] );
        error_tag_data[n] = target_tag_data[n] - gold_value;
        error_l2_norm += error_tag_data[n] * error_tag_data[n];
        tag_l2_norm += target_tag_data[n] * target_tag_data[n];
        error_tag_data[n] /= gold_value;
    }

    error_l2_norm = std::sqrt( error_l2_norm );
    tag_l2_norm = std::sqrt( tag_l2_norm );
    std::cout << "|e|_2 / |f|_2: " << error_l2_norm / tag_l2_norm << std::endl;

    error = target_iface->tag_set_data(
        target_error_tag, target_nodes.data(), num_target_nodes,
        static_cast<void *>( error_tag_data.getRawPtr() ) );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );

    // SOURCE MESH WRITE
    // -----------------

    error = source_iface->write_file( source_mesh_output_file.c_str(), 0, 0,
                                      &source_set, 1, &source_data_tag, 1 );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );

    // TARGET MESH WRITE
    // -----------------

    Teuchos::Array<moab::Tag> out_tags( 2 );
    out_tags[0] = target_data_tag;
    out_tags[1] = target_error_tag;
    error = target_iface->write_file( target_mesh_output_file.c_str(), 0, 0,
                                      &target_set, 1, &out_tags[0], 2 );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );
}

//---------------------------------------------------------------------------//
// end tstMoab_Mesh.cpp
//---------------------------------------------------------------------------//
