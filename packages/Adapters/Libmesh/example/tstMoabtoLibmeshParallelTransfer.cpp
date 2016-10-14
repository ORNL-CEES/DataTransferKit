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
 * \file   tstMoabTagFieldParallelPartition.cpp
 * \author Alex McCaskey
 * \brief  Tag field test
 */
//---------------------------------------------------------------------------//

#include <moab/Core.hpp>
#include <moab/Interface.hpp>
#include <moab/ParallelComm.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include "DTK_LibmeshManager.hpp"
#include "DTK_MapOperatorFactory.hpp"
#include "DTK_MoabManager.hpp"

#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_VerboseObject.hpp>

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

    // Extract the raw mpi communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<const Teuchos::MpiComm<int>> mpi_comm =
        Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>( comm );
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> opaque_comm =
        mpi_comm->getRawMpiComm();
    MPI_Comm raw_comm = ( *opaque_comm )();

    // Build the parameter list from the xml input.
    Teuchos::RCP<Teuchos::ParameterList> plist =
        Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( "input.xml",
                                          Teuchos::inoutArg( *plist ) );

    // Load the mesh.
    Teuchos::RCP<moab::Interface> source_iface =
        Teuchos::rcp( new moab::Core() );
    std::string meshFileName( "sahex1_unic.h5m" );
    std::string options =
        std::string( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;"
                     "PARALLEL_RESOLVE_SHARED_ENTS;PARTITION_DISTRIBUTE;"
                     "PARALLEL_GHOSTS=3.0.1;" );
    source_iface->load_file( meshFileName.c_str(), 0, options.c_str() );

    // Get the parallel moab instance.
    Teuchos::RCP<moab::ParallelComm> source_mesh = Teuchos::rcp(
        moab::ParallelComm::get_pcomm( source_iface.getRawPtr(), 0 ), false );

    // Get the entity set for the source part. Just use the root set for now.
    moab::EntityHandle source_set = source_iface->get_root_set();

    // Get the nodes in the source set.
    std::vector<moab::EntityHandle> source_nodes;
    moab::ErrorCode error = source_iface->get_entities_by_type(
        source_set, moab::MBVERTEX, source_nodes );
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

    // Make the tag parallel consistent.
    moab::Range shared_entities;
    error = source_mesh->get_shared_entities( -1, shared_entities, -1, false,
                                              true );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );
    error = source_mesh->exchange_tags( source_data_tag, shared_entities );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );

    // TARGET MESH READ
    // ----------------

    // Initialize libmesh.
    libMesh::LibMeshInit libmesh_init( argc, argv, raw_comm );

    // Create a target mesh.
    Teuchos::RCP<libMesh::Mesh> tgt_mesh =
        Teuchos::rcp( new libMesh::Mesh( libmesh_init.comm() ) );
    tgt_mesh->read( "sahex.e" );

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
    tgt_equation_systems.init();

    // SETUP DTK DATA
    // --------------

    // Create a manager for the source set elements.
    DataTransferKit::MoabManager src_manager( source_mesh, source_set );

    // Create a manager for the target set nodes.
    DataTransferKit::LibmeshManager tgt_manager(
        tgt_mesh, Teuchos::rcpFromRef( tgt_system ) );

    // Create a solution vector for the source.
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        src_vector = src_manager.createFieldMultiVector( source_node_set,
                                                         source_data_tag );

    // Create a solution vector for the target.
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        tgt_vector = tgt_manager.createFieldMultiVector( tgt_var_name );

    // Print out mesh info.
    Teuchos::RCP<Teuchos::FancyOStream> fancy_out =
        Teuchos::VerboseObjectBase::getDefaultOStream();
    fancy_out->setShowProcRank( true );
    src_manager.functionSpace()->entitySet()->describe( *fancy_out );
    tgt_manager.functionSpace()->entitySet()->describe( *fancy_out );

    // Create a map operator.
    Teuchos::ParameterList &dtk_list = plist->sublist( "DataTransferKit" );
    DataTransferKit::MapOperatorFactory op_factory;
    Teuchos::RCP<DataTransferKit::MapOperator> map_op = op_factory.create(
        src_vector->getMap(), tgt_vector->getMap(), dtk_list );

    // Setup the map operator.
    map_op->setup( src_manager.functionSpace(), tgt_manager.functionSpace() );

    // SOLUTION TRANSFER
    // -----------------

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
    libMesh::Mesh::const_element_iterator elements_begin =
        tgt_mesh->local_elements_begin();
    libMesh::Mesh::const_element_iterator elements_end =
        tgt_mesh->local_elements_end();
    std::vector<libMesh::dof_id_type> elem_dof_ids;
    int elem_n_nodes = 0;
    libMesh::Node *node;
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
            if ( comm->getRank() == node->processor_id() )
            {
                gold_value = dataFunction( ( *node )( 0 ), ( *node )( 1 ),
                                           ( *node )( 2 ) );
                tgt_value = tgt_system.solution->el( elem_dof_ids[n] );
                tgt_system.solution->set(
                    err_dof_ids[n], ( tgt_value - gold_value ) / gold_value );
                tgt_val_map.emplace( err_dof_ids[n], tgt_value * tgt_value );
                error_val_map.emplace( err_dof_ids[n],
                                       ( tgt_value - gold_value ) *
                                           ( tgt_value - gold_value ) );
            }
        }
    }
    tgt_system.solution->close();
    tgt_system.update();

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
    error = source_iface->write_file( "source_moab_out.h5m", "H5M",
                                      "PARALLEL=WRITE_PART", &source_set, 1,
                                      &source_data_tag, 1 );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );

    // TARGET MESH WRITE
    // -----------------

    libMesh::ExodusII_IO( *tgt_mesh )
        .write_equation_systems( "target_libmesh_out.exo",
                                 tgt_equation_systems );
}
