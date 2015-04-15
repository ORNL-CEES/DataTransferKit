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
 * \file   tstSTKMeshNodalShapeFunction.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  Nodal shape function test.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "DTK_STKMeshNodalShapeFunction.hpp"
#include "DTK_STKMeshEntity.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include <Teuchos_DefaultMpiComm.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_topology/topology.hpp>

//---------------------------------------------------------------------------//
// Hex-8 test.
TEUCHOS_UNIT_TEST( STKMeshNodalShapeFunction, hex_8_test )
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

    // Create a DTK entity for the hex.
    DataTransferKit::Entity dtk_entity = 
	DataTransferKit::STKMeshEntity( hex_entity, bulk_data.ptr() );

    // Create a shape function.
    Teuchos::RCP<DataTransferKit::EntityShapeFunction> shape_function =
	Teuchos::rcp( new DataTransferKit::STKMeshNodalShapeFunction(bulk_data) );

    // Test the shape function dof ids for the hex.
    Teuchos::Array<std::size_t> dof_ids;
    shape_function->entitySupportIds( dtk_entity, dof_ids );
    TEST_EQUALITY( num_nodes, dof_ids.size() );
    for ( unsigned n = 0; n < num_nodes; ++n )
    {
	TEST_EQUALITY( dof_ids[n], num_nodes*comm_rank + n + 5 );
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
	DataTransferKit::Entity dtk_node =
	    DataTransferKit::STKMeshEntity( nodes[n], bulk_data.ptr() );
	shape_function->entitySupportIds( dtk_node, dof_ids );
	TEST_EQUALITY( dof_ids.size(), 1 );
	TEST_EQUALITY( dof_ids[0], node_ids[n] );
    }
}

//---------------------------------------------------------------------------//
// end of tstSTKMeshNodalShapeFunction.cpp
//---------------------------------------------------------------------------//
