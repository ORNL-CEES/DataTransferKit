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
 * \file   tstSTKMeshEntityIntegrationRule.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  Integration rule function test.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "DTK_STKMeshEntityIntegrationRule.hpp"
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
TEUCHOS_UNIT_TEST( STKMeshEntityIntegrationRule, hex_8_test )
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

    // Create an integration rule.
    Teuchos::RCP<DataTransferKit::EntityIntegrationRule> integration_rule =
        Teuchos::rcp( new DataTransferKit::STKMeshEntityIntegrationRule(bulk_data) );

    // Test the integration rule.
    Teuchos::Array<Teuchos::Array<double> > p_1;
    Teuchos::Array<double> w_1;
    integration_rule->getIntegrationRule( dtk_entity, 1, p_1, w_1 );
    TEST_EQUALITY( 1, w_1.size() );
    TEST_EQUALITY( 1, p_1.size() );
    TEST_EQUALITY( 3, p_1[0].size() );
    TEST_EQUALITY( 8.0, w_1[0] );
    TEST_EQUALITY( 0.0, p_1[0][0] );
    TEST_EQUALITY( 0.0, p_1[0][1] );
    TEST_EQUALITY( 0.0, p_1[0][2] );

    Teuchos::Array<Teuchos::Array<double> > p_2;
    Teuchos::Array<double> w_2;
    integration_rule->getIntegrationRule( dtk_entity, 2, p_2, w_2 );
    TEST_EQUALITY( 8, w_2.size() );
    TEST_EQUALITY( 8, p_2.size() );
    for ( int i = 0; i < 8; ++i )
    {
        TEST_EQUALITY( w_2[i], 1.0 );
        TEST_EQUALITY( p_2[i].size(), 3 );

        for ( int d = 0; d < 3; ++d )
        {
            TEST_FLOATING_EQUALITY(
                std::abs(p_2[i][d]), 1.0 / std::sqrt(3.0), 1.0e-15 );
        }
    }
}

//---------------------------------------------------------------------------//
// end of tstSTKMeshEntityIntegrationRule.cpp
//---------------------------------------------------------------------------//
