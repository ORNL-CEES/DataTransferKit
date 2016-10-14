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
 * \file tstSTKMeshEntityLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief STKMeshEntityLocalMap unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_STKMeshEntity.hpp>
#include <DTK_STKMeshEntityLocalMap.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template <class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal>> getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp( new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Hex-8 test.
TEUCHOS_UNIT_TEST( STKMeshEntity, hex_8_test )
{
    // Extract the raw mpi communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm<int>();
    Teuchos::RCP<const Teuchos::MpiComm<int>> mpi_comm =
        Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>( comm );
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm>> opaque_comm =
        mpi_comm->getRawMpiComm();
    MPI_Comm raw_comm = ( *opaque_comm )();

    // Create meta data.
    int space_dim = 3;
    stk::mesh::MetaData meta_data( space_dim );

    // Make two parts.
    std::string p1_name = "part_1";
    stk::mesh::Part &part_1 = meta_data.declare_part( p1_name );
    stk::mesh::set_topology( part_1, stk::topology::HEX_8 );

    // Make a coordinate field.
    stk::mesh::Field<double, stk::mesh::Cartesian3d> &coord_field =
        meta_data
            .declare_field<stk::mesh::Field<double, stk::mesh::Cartesian3d>>(
                stk::topology::NODE_RANK, "coordinates" );
    meta_data.set_coordinate_field( &coord_field );
    stk::mesh::put_field( coord_field, part_1 );
    meta_data.commit();

    // Create bulk data.
    Teuchos::RCP<stk::mesh::BulkData> bulk_data =
        Teuchos::rcp( new stk::mesh::BulkData( meta_data, raw_comm ) );
    bulk_data->modification_begin();

    // Make a hex-8.
    int comm_rank = comm->getRank();
    stk::mesh::EntityId hex_id = 23 + comm_rank;
    stk::mesh::Entity hex_entity =
        bulk_data->declare_entity( stk::topology::ELEM_RANK, hex_id, part_1 );
    int num_nodes = 8;
    Teuchos::Array<stk::mesh::EntityId> node_ids( num_nodes );
    Teuchos::Array<stk::mesh::Entity> nodes( num_nodes );
    for ( int i = 0; i < num_nodes; ++i )
    {
        node_ids[i] = num_nodes * comm_rank + i + 5;
        nodes[i] = bulk_data->declare_entity( stk::topology::NODE_RANK,
                                              node_ids[i], part_1 );
        bulk_data->declare_relation( hex_entity, nodes[i], i );
    }
    bulk_data->modification_end();

    // Create the node coordinates.
    double *node_coords = 0;
    node_coords = stk::mesh::field_data( coord_field, nodes[0] );
    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[1] );
    node_coords[0] = 2.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[2] );
    node_coords[0] = 2.0;
    node_coords[1] = 2.0;
    node_coords[2] = 0.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[3] );
    node_coords[0] = 0.0;
    node_coords[1] = 2.0;
    node_coords[2] = 0.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[4] );
    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = 2.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[5] );
    node_coords[0] = 2.0;
    node_coords[1] = 0.0;
    node_coords[2] = 2.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[6] );
    node_coords[0] = 2.0;
    node_coords[1] = 2.0;
    node_coords[2] = 2.0;

    node_coords = stk::mesh::field_data( coord_field, nodes[7] );
    node_coords[0] = 0.0;
    node_coords[1] = 2.0;
    node_coords[2] = 2.0;

    // Create a local map from the bulk data.
    Teuchos::RCP<DataTransferKit::EntityLocalMap> local_map =
        Teuchos::rcp( new DataTransferKit::STKMeshEntityLocalMap( bulk_data ) );

    // Create a DTK entity for the hex.
    DataTransferKit::Entity dtk_entity =
        DataTransferKit::STKMeshEntity( hex_entity, bulk_data.ptr() );

    // Test the measure.
    TEST_EQUALITY( local_map->measure( dtk_entity ), 8.0 );

    // Test the centroid.
    Teuchos::Array<double> centroid( space_dim, 0.0 );
    local_map->centroid( dtk_entity, centroid() );
    TEST_EQUALITY( centroid[0], 1.0 );
    TEST_EQUALITY( centroid[1], 1.0 );
    TEST_EQUALITY( centroid[2], 1.0 );

    // Make a good point and a bad point.
    Teuchos::Array<double> good_point( space_dim );
    good_point[0] = 0.5;
    good_point[1] = 1.5;
    good_point[2] = 1.0;
    Teuchos::Array<double> bad_point( space_dim );
    bad_point[0] = 0.75;
    bad_point[1] = -1.75;
    bad_point[2] = 0.35;

    // Test the reference frame safeguard.
    TEST_ASSERT(
        local_map->isSafeToMapToReferenceFrame( dtk_entity, good_point() ) );
    TEST_ASSERT(
        !local_map->isSafeToMapToReferenceFrame( dtk_entity, bad_point() ) );

    // Test the mapping to reference frame.
    Teuchos::Array<double> ref_good_point( space_dim );
    bool good_map = local_map->mapToReferenceFrame( dtk_entity, good_point(),
                                                    ref_good_point() );
    TEST_ASSERT( good_map );
    TEST_EQUALITY( ref_good_point[0], -0.5 );
    TEST_EQUALITY( ref_good_point[1], 0.5 );
    TEST_EQUALITY( ref_good_point[2], 0.0 );

    Teuchos::Array<double> ref_bad_point( space_dim );
    bool bad_map = local_map->mapToReferenceFrame( dtk_entity, bad_point(),
                                                   ref_bad_point() );
    TEST_ASSERT( bad_map );

    // Test the point inclusion.
    TEST_ASSERT(
        local_map->checkPointInclusion( dtk_entity, ref_good_point() ) );
    TEST_ASSERT(
        !local_map->checkPointInclusion( dtk_entity, ref_bad_point() ) );

    // Test the map to physical frame.
    Teuchos::Array<double> phy_good_point( space_dim );
    local_map->mapToPhysicalFrame( dtk_entity, ref_good_point(),
                                   phy_good_point() );
    TEST_EQUALITY( good_point[0], phy_good_point[0] );
    TEST_EQUALITY( good_point[1], phy_good_point[1] );
    TEST_EQUALITY( good_point[2], phy_good_point[2] );

    Teuchos::Array<double> phy_bad_point( space_dim );
    local_map->mapToPhysicalFrame( dtk_entity, ref_bad_point(),
                                   phy_bad_point() );
    TEST_EQUALITY( bad_point[0], phy_bad_point[0] );
    TEST_EQUALITY( bad_point[1], phy_bad_point[1] );
    TEST_EQUALITY( bad_point[2], phy_bad_point[2] );

    // Test the coordinates of the points extracted through the centroid
    // function.
    DataTransferKit::Entity dtk_node;
    Teuchos::Array<double> point_coords( space_dim );
    for ( int n = 0; n < num_nodes; ++n )
    {
        dtk_node = DataTransferKit::STKMeshEntity( nodes[n], bulk_data.ptr() );
        local_map->centroid( dtk_node, point_coords() );
        node_coords = stk::mesh::field_data( coord_field, nodes[n] );
        TEST_EQUALITY( node_coords[0], point_coords[0] );
        TEST_EQUALITY( node_coords[1], point_coords[1] );
        TEST_EQUALITY( node_coords[2], point_coords[2] );
    }
}

//---------------------------------------------------------------------------//
// end tstSTKMeshEntityLocalMap.cpp
//---------------------------------------------------------------------------//
