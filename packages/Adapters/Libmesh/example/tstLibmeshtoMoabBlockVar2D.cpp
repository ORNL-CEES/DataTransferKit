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

    // Initialize libmesh.
    libMesh::LibMeshInit libmesh_init( argc, argv, raw_comm );

    // Create a target mesh.
    Teuchos::RCP<libMesh::Mesh> src_mesh =
        Teuchos::rcp( new libMesh::Mesh( libmesh_init.comm() ) );
    src_mesh->read( "proteus_fumex_2d.exo" );

    // Make a libmesh system. We will put a first order linear basis on the
    // elements for all subdomains.
    std::string src_var_name = "src_var";
    std::string src_err_name = "src_err";
    libMesh::EquationSystems src_equation_systems( *src_mesh );
    libMesh::System &src_system =
        src_equation_systems.add_system<libMesh::ExplicitSystem>(
            "Source System" );
    int src_var_id = src_system.add_variable( src_var_name );
    src_equation_systems.init();

    // Create some data and put it on the variable.
    libMesh::Mesh::const_element_iterator elements_begin =
        src_mesh->local_elements_begin();
    libMesh::Mesh::const_element_iterator elements_end =
        src_mesh->local_elements_end();
    int elem_n_nodes = 0;
    libMesh::Node *node;
    std::vector<libMesh::dof_id_type> elem_dof_ids;
    for ( auto element = elements_begin; element != elements_end; ++element )
    {
        src_system.get_dof_map().dof_indices( *element, elem_dof_ids,
                                              src_var_id );
        elem_n_nodes = ( *element )->n_nodes();

        double solution_val = ( *element )->subdomain_id() == 1 ? 10.0 : 2.0;

        for ( int n = 0; n < elem_n_nodes; ++n )
        {
            node = ( *element )->get_node( n );

            src_system.solution->set( elem_dof_ids[n], solution_val );
            //                dataFunction( (*node)(0), (*node)(1), (*node)(2))
            //                );
        }
    }
    src_system.solution->close();
    src_system.update();

    // Load the moab mesh.
    Teuchos::RCP<moab::Interface> target_iface =
        Teuchos::rcp( new moab::Core() );
    std::string meshFileName( "proteus_fumex_2d.h5m" );
    std::string options =
        std::string( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;"
                     "PARALLEL_RESOLVE_SHARED_ENTS;PARTITION_DISTRIBUTE;"
                     "PARALLEL_GHOSTS=3.0.1;" );
    target_iface->load_file( meshFileName.c_str(), 0, options.c_str() );
    moab::ErrorCode error = target_iface->set_dimension( 2 );
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
    double default_val = 0.0;
    bool created = false;
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

    // Create a manager for the source.
    DataTransferKit::LibmeshManager src_manager(
        src_mesh, Teuchos::rcpFromRef( src_system ) );

    // Create a manager for the target.
    DataTransferKit::MoabManager tgt_manager( target_mesh, target_set );

    // Create a solution vector for the source.
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        src_vector = src_manager.createFieldMultiVector( src_var_name );

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
    /*
        double gold_value = 0.0;
        double error_l2_norm = 0.0;
        double tag_l2_norm = 0.0;
        Teuchos::Array<double> error_tag_data( num_target_nodes );
        Teuchos::Array<double> target_coords( 3 );

        error = target_iface->tag_get_data( target_data_tag,
                                            target_nodes.data(),
                                            num_target_nodes,
                                            static_cast<void*>(target_tag_data.getRawPtr())
       );
        checkMoabErrorCode( error );
        assert( moab::MB_SUCCESS == error );

        for ( int n = 0; n < num_target_nodes; ++n )
        {
            error = target_iface->get_coords( &target_nodes[n],
                                              1,
                                              target_coords.getRawPtr() );
            checkMoabErrorCode( error );
            assert( moab::MB_SUCCESS == error );
             p
            double gold_value (*element->subdomain_id() == 1 ? 10.0 : 2.0);
            gold_value = dataFunction( target_coords[0],
                                       target_coords[1],
                                       target_coords[2] );
            error_tag_data[n] = target_tag_data[n] - gold_value;
            error_l2_norm += error_tag_data[n] * error_tag_data[n];
            tag_l2_norm += target_tag_data[n] * target_tag_data[n];
            error_tag_data[n] /= gold_value;
        }

        error_l2_norm = std::sqrt( error_l2_norm );
        tag_l2_norm = std::sqrt( tag_l2_norm );
        std::cout << "|e|_2 / |f|_2: " << error_l2_norm / tag_l2_norm <<
       std::endl;

        error = target_iface->tag_set_data( target_error_tag,
                                            target_nodes.data(),
                                            num_target_nodes,
                                            static_cast<void*>(error_tag_data.getRawPtr())
       );
        checkMoabErrorCode( error );
        assert( moab::MB_SUCCESS == error );

    */
    // SOURCE MESH WRITE
    // -----------------

    libMesh::ExodusII_IO( *src_mesh )
        .write_equation_systems( "source_libmesh_out.exo",
                                 src_equation_systems );

    // TARGET MESH WRITE
    // -----------------
    Teuchos::Array<moab::Tag> out_tags( 2 );
    out_tags[0] = target_data_tag;
    out_tags[1] = target_error_tag;
    error = target_iface->write_file( "target_moab_out.vtk", "VTK",
                                      "PARALLEL=WRITE_PART", &target_set, 1,
                                      &out_tags[0], 2 );
    checkMoabErrorCode( error );
    assert( moab::MB_SUCCESS == error );
}
