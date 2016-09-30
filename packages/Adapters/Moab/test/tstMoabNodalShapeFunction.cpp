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
 * \file   tstMoabNodalShapeFunction.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  Nodal shape function test.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "DTK_MoabNodalShapeFunction.hpp"
#include "DTK_MoabEntity.hpp"
#include <DTK_MoabHelpers.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include <Teuchos_DefaultMpiComm.hpp>

#include <moab/Interface.hpp>
#include <moab/ParallelComm.hpp>
#include <moab/Core.hpp>

//---------------------------------------------------------------------------//
// Hex-8 test.
TEUCHOS_UNIT_TEST( MoabNodalShapeFunction, hex_8_test )
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
    int space_dim = 3;
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

    node_coords[0] = 2.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;
    error = moab_mesh->create_vertex( node_coords, nodes[1] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 2.0;
    node_coords[1] = 2.0;
    node_coords[2] = 0.0;
    error = moab_mesh->create_vertex( node_coords, nodes[2] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 0.0;
    node_coords[1] = 2.0;
    node_coords[2] = 0.0;
    error = moab_mesh->create_vertex( node_coords, nodes[3] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = 2.0;
    error = moab_mesh->create_vertex( node_coords, nodes[4] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 2.0;
    node_coords[1] = 0.0;
    node_coords[2] = 2.0;
    error = moab_mesh->create_vertex( node_coords, nodes[5] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 2.0;
    node_coords[1] = 2.0;
    node_coords[2] = 2.0;
    error = moab_mesh->create_vertex( node_coords, nodes[6] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 0.0;
    node_coords[1] = 2.0;
    node_coords[2] = 2.0;
    error = moab_mesh->create_vertex( node_coords, nodes[7] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    // Make a hex-8.
    moab::EntityHandle hex_entity;
    error = moab_mesh->create_element( moab::MBHEX,
                                       nodes.getRawPtr(),
                                       8,
                                       hex_entity );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    // Index the sets in the mesh.
    Teuchos::RCP<DataTransferKit::MoabMeshSetIndexer> set_indexer =
        Teuchos::rcp( new DataTransferKit::MoabMeshSetIndexer(parallel_mesh) );

    // Create a DTK entity for the hex.
    DataTransferKit::Entity dtk_entity = DataTransferKit::MoabEntity(
        hex_entity, parallel_mesh.ptr(), set_indexer.ptr() );

    // Create a shape function.
    Teuchos::RCP<DataTransferKit::EntityShapeFunction> shape_function =
        Teuchos::rcp( new DataTransferKit::MoabNodalShapeFunction(parallel_mesh) );

    // Test the shape function dof ids for the hex.
    std::vector<DataTransferKit::EntityId> node_ids( num_nodes );
    DataTransferKit::MoabHelpers::getGlobalIds( *parallel_mesh,
                                                nodes.getRawPtr(),
                                                num_nodes,
                                                node_ids.data() );
    Teuchos::Array<DataTransferKit::SupportId> dof_ids;
    shape_function->entitySupportIds( dtk_entity, dof_ids );
    TEST_EQUALITY( num_nodes, dof_ids.size() );
    for ( unsigned n = 0; n < num_nodes; ++n )
    {
        TEST_EQUALITY( dof_ids[n], node_ids[n] );
    }

    // Test the value evaluation for the hex.
    Teuchos::Array<double> ref_point( space_dim, 0.0 );
    Teuchos::Array<double> values;
    shape_function->evaluateValue( dtk_entity, ref_point(), values );
    TEST_EQUALITY( values.size(), num_nodes );
    for ( unsigned n = 0; n < num_nodes; ++n )
    {
        TEST_EQUALITY( values[n], 1.0 / num_nodes );
    }
    ref_point[0] = -1.0;
    ref_point[1] = -1.0;
    ref_point[2] = -1.0;
    shape_function->evaluateValue( dtk_entity, ref_point(), values );
    TEST_EQUALITY( values.size(), num_nodes );
    TEST_EQUALITY( values[0], 1.0 );
    for ( unsigned n = 1; n < num_nodes; ++n )
    {
        TEST_EQUALITY( values[n], 0.0 );
    }

    // Test the gradient evaluation for the hex.
    Teuchos::Array<Teuchos::Array<double> > grads;
    ref_point.assign( 3, 0.0 );
    shape_function->evaluateGradient( dtk_entity, ref_point(), grads );
    TEST_EQUALITY( grads.size(), num_nodes );
    for ( unsigned n = 0; n < num_nodes; ++n )
    {
        TEST_EQUALITY( Teuchos::as<int>(grads[n].size()), space_dim );
    }

    TEST_EQUALITY( grads[0][0], -1.0 / num_nodes );
    TEST_EQUALITY( grads[0][1], -1.0 / num_nodes );
    TEST_EQUALITY( grads[0][2], -1.0 / num_nodes );

    TEST_EQUALITY( grads[1][0], 1.0 / num_nodes );
    TEST_EQUALITY( grads[1][1], -1.0 / num_nodes );
    TEST_EQUALITY( grads[1][2], -1.0 / num_nodes );

    TEST_EQUALITY( grads[2][0], 1.0 / num_nodes );
    TEST_EQUALITY( grads[2][1], 1.0 / num_nodes );
    TEST_EQUALITY( grads[2][2], -1.0 / num_nodes );

    TEST_EQUALITY( grads[3][0], -1.0 / num_nodes );
    TEST_EQUALITY( grads[3][1], 1.0 / num_nodes );
    TEST_EQUALITY( grads[3][2], -1.0 / num_nodes );

    TEST_EQUALITY( grads[4][0], -1.0 / num_nodes );
    TEST_EQUALITY( grads[4][1], -1.0 / num_nodes );
    TEST_EQUALITY( grads[4][2], 1.0 / num_nodes );

    TEST_EQUALITY( grads[5][0], 1.0 / num_nodes );
    TEST_EQUALITY( grads[5][1], -1.0 / num_nodes );
    TEST_EQUALITY( grads[5][2], 1.0 / num_nodes );

    TEST_EQUALITY( grads[6][0], 1.0 / num_nodes );
    TEST_EQUALITY( grads[6][1], 1.0 / num_nodes );
    TEST_EQUALITY( grads[6][2], 1.0 / num_nodes );

    TEST_EQUALITY( grads[7][0], -1.0 / num_nodes );
    TEST_EQUALITY( grads[7][1], 1.0 / num_nodes );
    TEST_EQUALITY( grads[7][2], 1.0 / num_nodes );

    // Test the shape function dof ids for the nodes.
    for ( unsigned n = 0; n < num_nodes; ++n )
    {
        dof_ids.clear();
        DataTransferKit::Entity dtk_node = DataTransferKit::MoabEntity(
            nodes[n], parallel_mesh.ptr(), set_indexer.ptr() );
        shape_function->entitySupportIds( dtk_node, dof_ids );
        TEST_EQUALITY( dof_ids.size(), 1 );
        TEST_EQUALITY( dof_ids[0], node_ids[n] );
    }
}

//---------------------------------------------------------------------------//
// Quad-4 test.
TEUCHOS_UNIT_TEST( MoabNodalShapeFunction, quad_4_test )
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
    int space_dim = 2;
    Teuchos::RCP<moab::Interface> moab_mesh = Teuchos::rcp( new moab::Core() );
    moab_mesh->set_dimension( space_dim );
    Teuchos::RCP<moab::ParallelComm> parallel_mesh =
        Teuchos::rcp( new moab::ParallelComm(moab_mesh.getRawPtr(),raw_comm) );

    // Create the nodes.
    moab::ErrorCode error = moab::MB_SUCCESS;
    unsigned num_nodes = 4;
    Teuchos::Array<moab::EntityHandle> nodes(num_nodes);
    double node_coords[2];
    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    error = moab_mesh->create_vertex( node_coords, nodes[0] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 2.0;
    node_coords[1] = 0.0;
    error = moab_mesh->create_vertex( node_coords, nodes[1] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 2.0;
    node_coords[1] = 2.0;
    error = moab_mesh->create_vertex( node_coords, nodes[2] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 0.0;
    node_coords[1] = 2.0;
    error = moab_mesh->create_vertex( node_coords, nodes[3] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    // Make a quad-4.
    moab::EntityHandle quad_entity;
    error = moab_mesh->create_element( moab::MBQUAD,
                                       nodes.getRawPtr(),
                                       4,
                                       quad_entity );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    // Index the sets in the mesh.
    Teuchos::RCP<DataTransferKit::MoabMeshSetIndexer> set_indexer =
        Teuchos::rcp( new DataTransferKit::MoabMeshSetIndexer(parallel_mesh) );

    // Create a DTK entity for the quad.
    DataTransferKit::Entity dtk_entity = DataTransferKit::MoabEntity(
        quad_entity, parallel_mesh.ptr(), set_indexer.ptr() );

    // Create a shape function.
    Teuchos::RCP<DataTransferKit::EntityShapeFunction> shape_function =
        Teuchos::rcp( new DataTransferKit::MoabNodalShapeFunction(parallel_mesh) );

    // Test the shape function dof ids for the quad.
    std::vector<DataTransferKit::EntityId> node_ids( num_nodes );
    DataTransferKit::MoabHelpers::getGlobalIds( *parallel_mesh,
                                                nodes.getRawPtr(),
                                                num_nodes,
                                                node_ids.data() );
    Teuchos::Array<DataTransferKit::SupportId> dof_ids;
    shape_function->entitySupportIds( dtk_entity, dof_ids );
    TEST_EQUALITY( num_nodes, dof_ids.size() );
    for ( unsigned n = 0; n < num_nodes; ++n )
    {
        TEST_EQUALITY( dof_ids[n], node_ids[n] );
    }

    // Test the value evaluation for the quad.
    Teuchos::Array<double> ref_point( space_dim, 0.0 );
    Teuchos::Array<double> values;
    shape_function->evaluateValue( dtk_entity, ref_point(), values );
    TEST_EQUALITY( values.size(), num_nodes );
    for ( unsigned n = 0; n < num_nodes; ++n )
    {
        TEST_EQUALITY( values[n], 1.0 / num_nodes );
    }
    ref_point[0] = -1.0;
    ref_point[1] = -1.0;
    shape_function->evaluateValue( dtk_entity, ref_point(), values );
    TEST_EQUALITY( values.size(), num_nodes );
    TEST_EQUALITY( values[0], 1.0 );
    for ( unsigned n = 1; n < num_nodes; ++n )
    {
        TEST_EQUALITY( values[n], 0.0 );
    }

    // Test the gradient evaluation for the quad.
    Teuchos::Array<Teuchos::Array<double> > grads;
    ref_point.assign( space_dim, 0.0 );
    shape_function->evaluateGradient( dtk_entity, ref_point(), grads );
    TEST_EQUALITY( grads.size(), num_nodes );
    for ( unsigned n = 0; n < num_nodes; ++n )
    {
        TEST_EQUALITY( Teuchos::as<int>(grads[n].size()), space_dim );
    }

    TEST_EQUALITY( grads[0][0], -1.0 / num_nodes );
    TEST_EQUALITY( grads[0][1], -1.0 / num_nodes );

    TEST_EQUALITY( grads[1][0], 1.0 / num_nodes );
    TEST_EQUALITY( grads[1][1], -1.0 / num_nodes );

    TEST_EQUALITY( grads[2][0], 1.0 / num_nodes );
    TEST_EQUALITY( grads[2][1], 1.0 / num_nodes );

    TEST_EQUALITY( grads[3][0], -1.0 / num_nodes );
    TEST_EQUALITY( grads[3][1], 1.0 / num_nodes );

    // Test the shape function dof ids for the nodes.
    for ( unsigned n = 0; n < num_nodes; ++n )
    {
        dof_ids.clear();
        DataTransferKit::Entity dtk_node = DataTransferKit::MoabEntity(
            nodes[n], parallel_mesh.ptr(), set_indexer.ptr() );
        shape_function->entitySupportIds( dtk_node, dof_ids );
        TEST_EQUALITY( dof_ids.size(), 1 );
        TEST_EQUALITY( dof_ids[0], node_ids[n] );
    }
}

//---------------------------------------------------------------------------//
// end of tstMoabNodalShapeFunction.cpp
//---------------------------------------------------------------------------//
