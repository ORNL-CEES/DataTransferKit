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
 * \file   tstLibmeshVariableField.cpp
 * \author Stuart Slattery
 * \brief  Libmesh variable vector test.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_FieldMultiVector.hpp>
#include <DTK_LibmeshManager.hpp>
#include <DTK_LibmeshVariableField.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include <libmesh/cell_hex8.h>
#include <libmesh/dof_map.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/libmesh.h>
#include <libmesh/linear_partitioner.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/node.h>
#include <libmesh/parallel.h>
#include <libmesh/point.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/system.h>

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( LibmeshVariableField, push_pull_test )
{
    // Extract the raw mpi communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    Teuchos::RCP<const Teuchos::MpiComm<int>> mpi_comm =
        Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>( comm );
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> opaque_comm =
        mpi_comm->getRawMpiComm();
    MPI_Comm raw_comm = ( *opaque_comm )();

    // Create the mesh.
    TEST_ASSERT( !libMesh::initialized() );
    const std::string argv_string = "unit_test";
    const char *argv_char = argv_string.c_str();
    libMesh::LibMeshInit libmesh_init( 1, &argv_char, raw_comm );
    TEST_ASSERT( libMesh::initialized() );
    TEST_EQUALITY( (int)libmesh_init.comm().rank(), comm->getRank() );
    Teuchos::RCP<libMesh::Mesh> mesh =
        Teuchos::rcp( new libMesh::Mesh( libmesh_init.comm() ) );
    libMesh::MeshTools::Generation::build_cube( *mesh, 4, 4, 4, 0.0, 1.0, 0.0,
                                                1.0, 0.0, 1.0, libMesh::HEX8 );

    // Parition the mesh.
    libMesh::LinearPartitioner partitioner;
    partitioner.partition( *mesh );

    // Check libmesh validity.
    mesh->libmesh_assert_valid_parallel_ids();

    // Make a libmesh system. We will put a first order linear basis on the
    // elements for all subdomains.
    std::string var_name = "test_var";
    libMesh::EquationSystems equation_systems( *mesh );
    libMesh::ExplicitSystem &system =
        equation_systems.add_system<libMesh::ExplicitSystem>( "Test System" );
    int var_id = system.add_variable( var_name );
    equation_systems.init();

    // Put some data in the variable.
    int sys_id = system.number();
    libMesh::MeshBase::node_iterator nodes_begin = mesh->local_nodes_begin();
    libMesh::MeshBase::node_iterator nodes_end = mesh->local_nodes_end();
    for ( auto node_it = nodes_begin; node_it != nodes_end; ++node_it )
    {
        auto dof_id = ( *node_it )->dof_number( sys_id, var_id, 0 );
        system.solution->set( dof_id, ( *node_it )->id() );
    }
    system.solution->close();
    system.update();

    // Create a vector from the variable over all subdomains.
    DataTransferKit::LibmeshManager manager( mesh,
                                             Teuchos::rcpFromRef( system ) );
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        var_vec = manager.createFieldMultiVector( var_name );

    // Test the vector size.
    unsigned int num_nodes = mesh->n_local_nodes();
    TEST_EQUALITY( 1, Teuchos::as<int>( var_vec->getNumVectors() ) );
    TEST_EQUALITY( num_nodes, var_vec->getLocalLength() );
    TEST_EQUALITY( 125, var_vec->getGlobalLength() );

    // Test the variable data and add some.
    Teuchos::rcp_dynamic_cast<DataTransferKit::FieldMultiVector>( var_vec )
        ->pullDataFromApplication();
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<double>> var_vec_view =
        var_vec->get2dViewNonConst();
    int n = 0;
    for ( auto node_it = nodes_begin; node_it != nodes_end; ++node_it )
    {
        TEST_EQUALITY( var_vec_view[0][n], ( *node_it )->id() );
        var_vec_view[0][n] += 1.0;
        ++n;
    }

    // Push the data back to Libmesh
    Teuchos::rcp_dynamic_cast<DataTransferKit::FieldMultiVector>( var_vec )
        ->pushDataToApplication();

    // Test the Libmesh variable.
    for ( auto node = nodes_begin; node != nodes_end; ++node )
    {
        auto dof_id = ( *node )->dof_number( sys_id, var_id, 0 );
        TEST_EQUALITY( system.solution->el( dof_id ), ( *node )->id() + 1 );
    }
}

//---------------------------------------------------------------------------//
// end of tstLibmeshVariableField.cpp
//---------------------------------------------------------------------------//
