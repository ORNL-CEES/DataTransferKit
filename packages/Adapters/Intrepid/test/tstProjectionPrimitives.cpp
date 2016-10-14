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
 * \file   tstProjectionPrimitives.cpp
 * \author Stuart Slattery
 * \brief  Projection primitives tests.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopology.hpp>

#include <Intrepid_FieldContainer.hpp>

#include "DTK_ProjectionPrimitives.hpp"

//---------------------------------------------------------------------------//
// DEFAULT TEST PARAMETERS.
//---------------------------------------------------------------------------//
Teuchos::ParameterList defaultParameters()
{
    Teuchos::ParameterList plist;
    plist.set<int>( "Max Newton Iterations", 1000 );
    plist.set<double>( "Newton Tolerance", 1.0e-9 );
    plist.set<double>( "Geometric Tolerance", 1.0e-6 );
    plist.set<double>( "Merge Tolerance", 1.0e-3 );
    plist.set<double>( "Spatial Tolerance", 1.0e-3 );
    plist.set<double>( "Normal Tolerance", std::sqrt( 2.0 ) / 2.0 );
    return plist;
}

//---------------------------------------------------------------------------//
// FLOATING POINT EPSILON FOR FLOATING EQUALITY CHECKS
//---------------------------------------------------------------------------//

const double epsilon = 1.0e-8;

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ProjectionPrimitives, quad_face_volume_of_influence_test )
{
    Teuchos::ParameterList parameters = defaultParameters();

    int elem_num_nodes = 4;
    unsigned space_dim = 3;
    Intrepid::FieldContainer<double> point( space_dim );
    Intrepid::FieldContainer<double> side_nodes( elem_num_nodes, space_dim );
    Intrepid::FieldContainer<double> side_node_normals( elem_num_nodes,
                                                        space_dim );
    shards::CellTopology side_topo =
        shards::getCellTopologyData<shards::Quadrilateral<4>>();
    bool in_volume_of_influence = false;

    side_nodes( 0, 0 ) = -1.0;
    side_nodes( 0, 1 ) = -1.0;
    side_nodes( 0, 2 ) = 0.0;
    side_nodes( 1, 0 ) = 1.0;
    side_nodes( 1, 1 ) = -1.0;
    side_nodes( 1, 2 ) = 0.0;
    side_nodes( 2, 0 ) = 1.0;
    side_nodes( 2, 1 ) = 1.0;
    side_nodes( 2, 2 ) = 0.0;
    side_nodes( 3, 0 ) = -1.0;
    side_nodes( 3, 1 ) = 1.0;
    side_nodes( 3, 2 ) = 0.0;

    side_node_normals( 0, 0 ) = 0.0;
    side_node_normals( 0, 1 ) = 0.0;
    side_node_normals( 0, 2 ) = 1.0;
    side_node_normals( 1, 0 ) = 0.0;
    side_node_normals( 1, 1 ) = 0.0;
    side_node_normals( 1, 2 ) = 1.0;
    side_node_normals( 2, 0 ) = 0.0;
    side_node_normals( 2, 1 ) = 0.0;
    side_node_normals( 2, 2 ) = 1.0;
    side_node_normals( 3, 0 ) = 0.0;
    side_node_normals( 3, 1 ) = 0.0;
    side_node_normals( 3, 2 ) = 1.0;

    // Inclusion test 1.
    point( 0 ) = 0.2;
    point( 1 ) = -0.34;
    point( 2 ) = 1.1;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( in_volume_of_influence );

    // Inclusion test 2.
    point( 0 ) = 2.2;
    point( 1 ) = -0.34;
    point( 2 ) = 1.1;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( !in_volume_of_influence );

    // Inclusion test 3.
    point( 0 ) = -1.0;
    point( 1 ) = -1.0;
    point( 2 ) = 0.0;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( in_volume_of_influence );

    // Inclusion test 4.
    point( 0 ) = -1.0;
    point( 1 ) = -1.0;
    point( 2 ) = 1.0;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( in_volume_of_influence );

    // Inclusion test 5.
    point( 0 ) = 1.0 + epsilon;
    point( 1 ) = 1.0 + epsilon;
    point( 2 ) = 0.0;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( in_volume_of_influence );

    // Inclusion test 6.
    point( 0 ) = 0.3;
    point( 1 ) = -1.0 - epsilon;
    point( 2 ) = 0.0;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( in_volume_of_influence );

    // Inclusion test 7.
    point( 0 ) = 0.333333;
    point( 1 ) = -1.0;
    point( 2 ) = 0.0;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( in_volume_of_influence );

    // Inclusion test 8.
    point( 0 ) = 2.3;
    point( 1 ) = 0.22;
    point( 2 ) = 0.0;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( !in_volume_of_influence );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ProjectionPrimitives, tri_face_volume_of_influence_test )
{
    Teuchos::ParameterList parameters = defaultParameters();

    int elem_num_nodes = 3;
    unsigned space_dim = 3;
    Intrepid::FieldContainer<double> point( space_dim );
    Intrepid::FieldContainer<double> side_nodes( elem_num_nodes, space_dim );
    Intrepid::FieldContainer<double> side_node_normals( elem_num_nodes,
                                                        space_dim );
    shards::CellTopology side_topo =
        shards::getCellTopologyData<shards::Triangle<3>>();
    bool in_volume_of_influence = false;

    side_nodes( 0, 0 ) = -1.0;
    side_nodes( 0, 1 ) = -1.0;
    side_nodes( 0, 2 ) = 0.0;
    side_nodes( 1, 0 ) = 1.0;
    side_nodes( 1, 1 ) = -1.0;
    side_nodes( 1, 2 ) = 0.0;
    side_nodes( 2, 0 ) = 1.0;
    side_nodes( 2, 1 ) = 1.0;
    side_nodes( 2, 2 ) = 0.0;

    side_node_normals( 0, 0 ) = 0.0;
    side_node_normals( 0, 1 ) = 0.0;
    side_node_normals( 0, 2 ) = 1.0;
    side_node_normals( 1, 0 ) = 0.0;
    side_node_normals( 1, 1 ) = 0.0;
    side_node_normals( 1, 2 ) = 1.0;
    side_node_normals( 2, 0 ) = 0.0;
    side_node_normals( 2, 1 ) = 0.0;
    side_node_normals( 2, 2 ) = 1.0;

    // Inclusion test 1.
    point( 0 ) = 0.2;
    point( 1 ) = -0.34;
    point( 2 ) = 1.1;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( in_volume_of_influence );

    // Inclusion test 2.
    point( 0 ) = 2.2;
    point( 1 ) = -0.34;
    point( 2 ) = 1.1;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( !in_volume_of_influence );

    // Inclusion test 3.
    point( 0 ) = -1.0;
    point( 1 ) = -1.0;
    point( 2 ) = 0.0;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( in_volume_of_influence );

    // Inclusion test 4.
    point( 0 ) = -1.0;
    point( 1 ) = -1.0;
    point( 2 ) = 1.0;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( in_volume_of_influence );

    // Inclusion test 5.
    point( 0 ) = 1.0 + epsilon;
    point( 1 ) = 1.0 + epsilon;
    point( 2 ) = 0.0;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( in_volume_of_influence );

    // Inclusion test 6.
    point( 0 ) = 0.3;
    point( 1 ) = -1.0 - epsilon;
    point( 2 ) = 0.0;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( in_volume_of_influence );

    // Inclusion test 7.
    point( 0 ) = 0.333333;
    point( 1 ) = -1.0;
    point( 2 ) = 0.0;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( in_volume_of_influence );

    // Inclusion test 8.
    point( 0 ) = 2.3;
    point( 1 ) = 0.22;
    point( 2 ) = 0.0;

    in_volume_of_influence =
        DataTransferKit::ProjectionPrimitives::pointInFaceVolumeOfInfluence(
            parameters, point, side_nodes, side_node_normals, side_topo );

    TEST_ASSERT( !in_volume_of_influence );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ProjectionPrimitives, quad_blue_to_green_projection_test )
{
    Teuchos::ParameterList parameters = defaultParameters();

    int elem_num_nodes = 4;
    unsigned space_dim = 3;
    Intrepid::FieldContainer<double> point( space_dim );
    Intrepid::FieldContainer<double> proj_point( space_dim );
    Intrepid::FieldContainer<double> param_point( 1, space_dim );
    Intrepid::FieldContainer<double> side_nodes( elem_num_nodes, space_dim );
    Intrepid::FieldContainer<double> side_node_normals( elem_num_nodes,
                                                        space_dim );
    shards::CellTopology side_topo =
        shards::getCellTopologyData<shards::Quadrilateral<4>>();
    int green_edge_id = 0;
    int green_node_id = 0;

    side_nodes( 0, 0 ) = -1.0;
    side_nodes( 0, 1 ) = -1.0;
    side_nodes( 0, 2 ) = 0.0;
    side_nodes( 1, 0 ) = 1.0;
    side_nodes( 1, 1 ) = -1.0;
    side_nodes( 1, 2 ) = 0.0;
    side_nodes( 2, 0 ) = 1.0;
    side_nodes( 2, 1 ) = 1.0;
    side_nodes( 2, 2 ) = 0.0;
    side_nodes( 3, 0 ) = -1.0;
    side_nodes( 3, 1 ) = 1.0;
    side_nodes( 3, 2 ) = 0.0;

    side_node_normals( 0, 0 ) = 0.0;
    side_node_normals( 0, 1 ) = 0.0;
    side_node_normals( 0, 2 ) = 1.0;
    side_node_normals( 1, 0 ) = 0.0;
    side_node_normals( 1, 1 ) = 0.0;
    side_node_normals( 1, 2 ) = 1.0;
    side_node_normals( 2, 0 ) = 0.0;
    side_node_normals( 2, 1 ) = 0.0;
    side_node_normals( 2, 2 ) = 1.0;
    side_node_normals( 3, 0 ) = 0.0;
    side_node_normals( 3, 1 ) = 0.0;
    side_node_normals( 3, 2 ) = 1.0;

    // Projection test 1.
    point( 0 ) = 0.2;
    point( 1 ) = -0.34;
    point( 2 ) = 1.1;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), 0.2 );
    TEST_FLOATING_EQUALITY( proj_point( 1 ), -0.34, epsilon );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, -1 );
    TEST_EQUALITY( green_node_id, -1 );

    // Projection test 2.
    point( 0 ) = -1.0;
    point( 1 ) = -1.0;
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), -1.0 );
    TEST_EQUALITY( proj_point( 1 ), -1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, -1 );
    TEST_EQUALITY( green_node_id, 0 );

    // Projection test 3.
    point( 0 ) = -1.0;
    point( 1 ) = -1.0;
    point( 2 ) = 1.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), -1.0 );
    TEST_EQUALITY( proj_point( 1 ), -1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, -1 );
    TEST_EQUALITY( green_node_id, 0 );

    // Projection test 4.
    point( 0 ) = 1.0 + epsilon;
    point( 1 ) = 1.0 + epsilon;
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), 1.0 );
    TEST_EQUALITY( proj_point( 1 ), 1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, -1 );
    TEST_EQUALITY( green_node_id, 2 );

    // Projection test 5.
    point( 0 ) = 0.3;
    point( 1 ) = -1.0 - epsilon;
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_FLOATING_EQUALITY( proj_point( 0 ), 0.3, epsilon );
    TEST_EQUALITY( proj_point( 1 ), -1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, 0 );
    TEST_EQUALITY( green_node_id, -1 );

    // Projection test 6.
    point( 0 ) = 0.333333;
    point( 1 ) = -1.0;
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_FLOATING_EQUALITY( proj_point( 0 ), 0.333333, epsilon );
    TEST_EQUALITY( proj_point( 1 ), -1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, 0 );
    TEST_EQUALITY( green_node_id, -1 );

    // Projection test 7.
    point( 0 ) = 1.0;
    point( 1 ) = -1.0;
    point( 2 ) = 1.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), 1.0 );
    TEST_EQUALITY( proj_point( 1 ), -1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, -1 );
    TEST_EQUALITY( green_node_id, 1 );

    // Projection test 8.
    point( 0 ) = 1.0;
    point( 1 ) = 1.0;
    point( 2 ) = 1.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), 1.0 );
    TEST_EQUALITY( proj_point( 1 ), 1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, -1 );
    TEST_EQUALITY( green_node_id, 2 );

    // Projection test 9.
    point( 0 ) = -1.0;
    point( 1 ) = 1.0;
    point( 2 ) = 1.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), -1.0 );
    TEST_EQUALITY( proj_point( 1 ), 1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, -1 );
    TEST_EQUALITY( green_node_id, 3 );

    // Projection test 10.
    point( 0 ) = 0.333333;
    point( 1 ) = 1.0;
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_FLOATING_EQUALITY( proj_point( 0 ), 0.333333, epsilon );
    TEST_EQUALITY( proj_point( 1 ), 1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, 2 );
    TEST_EQUALITY( green_node_id, -1 );

    // Projection test 11.
    point( 0 ) = 1.0;
    point( 1 ) = 0.333333;
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), 1.0 );
    TEST_FLOATING_EQUALITY( proj_point( 1 ), 0.333333, epsilon );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, 1 );
    TEST_EQUALITY( green_node_id, -1 );

    // Projection test 12.
    point( 0 ) = -1.0;
    point( 1 ) = 0.333333;
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), -1.0 );
    TEST_FLOATING_EQUALITY( proj_point( 1 ), 0.333333, epsilon );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, 3 );
    TEST_EQUALITY( green_node_id, -1 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ProjectionPrimitives, tri_blue_to_green_projection_test )
{
    int elem_num_nodes = 3;
    unsigned space_dim = 3;
    Intrepid::FieldContainer<double> point( space_dim );
    Intrepid::FieldContainer<double> param_point( 1, space_dim );
    Intrepid::FieldContainer<double> proj_point( space_dim );
    Intrepid::FieldContainer<double> side_nodes( elem_num_nodes, space_dim );
    Intrepid::FieldContainer<double> side_node_normals( elem_num_nodes,
                                                        space_dim );
    shards::CellTopology side_topo =
        shards::getCellTopologyData<shards::Triangle<3>>();
    int green_edge_id = 0;
    int green_node_id = 0;
    Teuchos::ParameterList parameters = defaultParameters();

    side_nodes( 0, 0 ) = -1.0;
    side_nodes( 0, 1 ) = -1.0;
    side_nodes( 0, 2 ) = 0.0;
    side_nodes( 1, 0 ) = 1.0;
    side_nodes( 1, 1 ) = -1.0;
    side_nodes( 1, 2 ) = 0.0;
    side_nodes( 2, 0 ) = 1.0;
    side_nodes( 2, 1 ) = 1.0;
    side_nodes( 2, 2 ) = 0.0;

    side_node_normals( 0, 0 ) = 0.0;
    side_node_normals( 0, 1 ) = 0.0;
    side_node_normals( 0, 2 ) = 1.0;
    side_node_normals( 1, 0 ) = 0.0;
    side_node_normals( 1, 1 ) = 0.0;
    side_node_normals( 1, 2 ) = 1.0;
    side_node_normals( 2, 0 ) = 0.0;
    side_node_normals( 2, 1 ) = 0.0;
    side_node_normals( 2, 2 ) = 1.0;

    // Projection test 1.
    point( 0 ) = 0.2;
    point( 1 ) = -0.34;
    point( 2 ) = 1.1;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_FLOATING_EQUALITY( proj_point( 0 ), 0.2,
                            parameters.get<double>( "Newton Tolerance" ) );
    TEST_FLOATING_EQUALITY( proj_point( 1 ), -0.34,
                            parameters.get<double>( "Newton Tolerance" ) );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, -1 );
    TEST_EQUALITY( green_node_id, -1 );

    // Projection test 2.
    point( 0 ) = -1.0;
    point( 1 ) = -1.0;
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), -1.0 );
    TEST_EQUALITY( proj_point( 1 ), -1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, -1 );
    TEST_EQUALITY( green_node_id, 0 );

    // Projection test 3.
    point( 0 ) = -1.0;
    point( 1 ) = -1.0;
    point( 2 ) = 1.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), -1.0 );
    TEST_EQUALITY( proj_point( 1 ), -1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, -1 );
    TEST_EQUALITY( green_node_id, 0 );

    // Projection test 4.
    point( 0 ) = 1.0 + parameters.get<double>( "Newton Tolerance" );
    point( 1 ) = 1.0 + parameters.get<double>( "Newton Tolerance" );
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), 1.0 );
    TEST_EQUALITY( proj_point( 1 ), 1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, -1 );
    TEST_EQUALITY( green_node_id, 2 );

    // Projection test 5.
    point( 0 ) = 0.3;
    point( 1 ) = -1.0 - parameters.get<double>( "Newton Tolerance" );
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_FLOATING_EQUALITY( proj_point( 0 ), 0.3,
                            10.0 *
                                parameters.get<double>( "Newton Tolerance" ) );
    TEST_EQUALITY( proj_point( 1 ), -1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, 0 );
    TEST_EQUALITY( green_node_id, -1 );

    // Projection test 6.
    point( 0 ) = 0.333333;
    point( 1 ) = -1.0;
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_FLOATING_EQUALITY( proj_point( 0 ), 0.333333, epsilon );
    TEST_EQUALITY( proj_point( 1 ), -1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, 0 );
    TEST_EQUALITY( green_node_id, -1 );

    // Projection test 7.
    point( 0 ) = 1.0 + parameters.get<double>( "Newton Tolerance" );
    point( 1 ) = -1.0 - parameters.get<double>( "Newton Tolerance" );
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), 1.0 );
    TEST_EQUALITY( proj_point( 1 ), -1.0 );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, -1 );
    TEST_EQUALITY( green_node_id, 1 );

    // Projection test 8.
    point( 0 ) = 1.0;
    point( 1 ) = 0.333333;
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_EQUALITY( proj_point( 0 ), 1.0 );
    TEST_FLOATING_EQUALITY( proj_point( 1 ), 0.333333, epsilon );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, 1 );
    TEST_EQUALITY( green_node_id, -1 );

    // Projection test 9.
    point( 0 ) = 0.2;
    point( 1 ) = 0.2;
    point( 2 ) = 0.0;

    DataTransferKit::ProjectionPrimitives::projectPointToFace(
        parameters, point, side_nodes, side_node_normals, side_topo,
        param_point, proj_point, green_edge_id, green_node_id );

    TEST_FLOATING_EQUALITY( proj_point( 0 ), 0.2, epsilon );
    TEST_FLOATING_EQUALITY( proj_point( 1 ), 0.2, epsilon );
    TEST_EQUALITY( proj_point( 2 ), 0.0 );
    TEST_EQUALITY( green_edge_id, 2 );
    TEST_EQUALITY( green_node_id, -1 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ProjectionPrimitives, project_blue_feature_test )
{
    typedef DataTransferKit::ProjectionPrimitives PP;

    Teuchos::ParameterList parameters = defaultParameters();

    // Allocate Arrays.
    int space_dim = 3;
    Intrepid::FieldContainer<double> blue_point( space_dim );
    Intrepid::FieldContainer<double> blue_normal( space_dim );
    Intrepid::FieldContainer<double> green_edge_nodes( 2, space_dim );
    Intrepid::FieldContainer<double> green_edge_node_normals( 2, space_dim );
    Intrepid::FieldContainer<double> projected_blue_point( space_dim );
    bool has_projection = false;
    int green_node_id = -1;

    // Projection 1.
    blue_point( 0 ) = 1.0;
    blue_point( 1 ) = 1.0;
    blue_point( 2 ) = 1.0;

    blue_normal( 0 ) = 0.0;
    blue_normal( 1 ) = 0.0;
    blue_normal( 2 ) = -1.0;

    green_edge_nodes( 0, 0 ) = -1.0;
    green_edge_nodes( 0, 1 ) = -1.0;
    green_edge_nodes( 0, 2 ) = 0.0;
    green_edge_nodes( 1, 0 ) = 2.0;
    green_edge_nodes( 1, 1 ) = 2.0;
    green_edge_nodes( 1, 2 ) = 0.0;

    green_edge_node_normals( 0, 0 ) = 0.0;
    green_edge_node_normals( 0, 1 ) = 0.0;
    green_edge_node_normals( 0, 2 ) = 1.0;
    green_edge_node_normals( 1, 0 ) = 0.0;
    green_edge_node_normals( 1, 1 ) = 0.0;
    green_edge_node_normals( 1, 2 ) = 1.0;

    has_projection = PP::projectPointFeatureToEdgeFeature(
        parameters, blue_point, blue_normal, green_edge_nodes,
        green_edge_node_normals, projected_blue_point, green_node_id );

    TEST_ASSERT( has_projection );
    TEST_EQUALITY( green_node_id, -1 );
    TEST_FLOATING_EQUALITY( projected_blue_point( 0 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( projected_blue_point( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( projected_blue_point( 2 ), 0.0, epsilon );

    // Projection 2.
    blue_point( 0 ) = 1.0;
    blue_point( 1 ) = 1.0;
    blue_point( 2 ) = 1.0;

    blue_normal( 0 ) = 0.0;
    blue_normal( 1 ) = 0.0;
    blue_normal( 2 ) = -1.0;

    green_edge_nodes( 0, 0 ) = 0.0;
    green_edge_nodes( 0, 1 ) = 0.0;
    green_edge_nodes( 0, 2 ) = 1.0;
    green_edge_nodes( 1, 0 ) = 2.0;
    green_edge_nodes( 1, 1 ) = 2.0;
    green_edge_nodes( 1, 2 ) = 1.0;

    green_edge_node_normals( 0, 0 ) = 0.0;
    green_edge_node_normals( 0, 1 ) = 0.0;
    green_edge_node_normals( 0, 2 ) = 1.0;
    green_edge_node_normals( 1, 0 ) = 0.0;
    green_edge_node_normals( 1, 1 ) = 0.0;
    green_edge_node_normals( 1, 2 ) = 1.0;

    has_projection = PP::projectPointFeatureToEdgeFeature(
        parameters, blue_point, blue_normal, green_edge_nodes,
        green_edge_node_normals, projected_blue_point, green_node_id );

    TEST_ASSERT( has_projection );
    TEST_EQUALITY( green_node_id, -1 );
    TEST_EQUALITY( projected_blue_point( 0 ), 1.0 );
    TEST_EQUALITY( projected_blue_point( 1 ), 1.0 );
    TEST_EQUALITY( projected_blue_point( 2 ), 1.0 );

    // Projection 3.
    blue_point( 0 ) = 2.0;
    blue_point( 1 ) = 2.0;
    blue_point( 2 ) = 1.0;

    blue_normal( 0 ) = 0.0;
    blue_normal( 1 ) = 0.0;
    blue_normal( 2 ) = -1.0;

    green_edge_nodes( 0, 0 ) = 0.0;
    green_edge_nodes( 0, 1 ) = 0.0;
    green_edge_nodes( 0, 2 ) = 1.0;
    green_edge_nodes( 1, 0 ) = 2.0;
    green_edge_nodes( 1, 1 ) = 2.0;
    green_edge_nodes( 1, 2 ) = 1.0;

    green_edge_node_normals( 0, 0 ) = 0.0;
    green_edge_node_normals( 0, 1 ) = 0.0;
    green_edge_node_normals( 0, 2 ) = 1.0;
    green_edge_node_normals( 1, 0 ) = 0.0;
    green_edge_node_normals( 1, 1 ) = 0.0;
    green_edge_node_normals( 1, 2 ) = 1.0;

    has_projection = PP::projectPointFeatureToEdgeFeature(
        parameters, blue_point, blue_normal, green_edge_nodes,
        green_edge_node_normals, projected_blue_point, green_node_id );

    TEST_ASSERT( has_projection );
    TEST_EQUALITY( green_node_id, 1 );
    TEST_EQUALITY( projected_blue_point( 0 ), 2.0 );
    TEST_EQUALITY( projected_blue_point( 1 ), 2.0 );
    TEST_EQUALITY( projected_blue_point( 2 ), 1.0 );

    // Projection 4.
    blue_point( 0 ) = 2.0;
    blue_point( 1 ) = 2.0;
    blue_point( 2 ) = 1.0;

    blue_normal( 0 ) = 0.0;
    blue_normal( 1 ) = 0.0;
    blue_normal( 2 ) = -1.0;

    green_edge_nodes( 0, 0 ) = 0.0;
    green_edge_nodes( 0, 1 ) = 0.0;
    green_edge_nodes( 0, 2 ) = 0.0;
    green_edge_nodes( 1, 0 ) = 2.0;
    green_edge_nodes( 1, 1 ) = 2.0;
    green_edge_nodes( 1, 2 ) = 0.0;

    green_edge_node_normals( 0, 0 ) = 0.0;
    green_edge_node_normals( 0, 1 ) = 0.0;
    green_edge_node_normals( 0, 2 ) = 1.0;
    green_edge_node_normals( 1, 0 ) = 0.0;
    green_edge_node_normals( 1, 1 ) = 0.0;
    green_edge_node_normals( 1, 2 ) = 1.0;

    has_projection = PP::projectPointFeatureToEdgeFeature(
        parameters, blue_point, blue_normal, green_edge_nodes,
        green_edge_node_normals, projected_blue_point, green_node_id );

    TEST_ASSERT( has_projection );
    TEST_EQUALITY( green_node_id, 1 );
    TEST_EQUALITY( projected_blue_point( 0 ), 2.0 );
    TEST_EQUALITY( projected_blue_point( 1 ), 2.0 );
    TEST_EQUALITY( projected_blue_point( 2 ), 0.0 );

    // Projection 5.
    blue_point( 0 ) = 0.0;
    blue_point( 1 ) = 0.0;
    blue_point( 2 ) = 1.0;

    blue_normal( 0 ) = 0.0;
    blue_normal( 1 ) = 0.0;
    blue_normal( 2 ) = -1.0;

    green_edge_nodes( 0, 0 ) = 0.0;
    green_edge_nodes( 0, 1 ) = 0.0;
    green_edge_nodes( 0, 2 ) = 0.0;
    green_edge_nodes( 1, 0 ) = 2.0;
    green_edge_nodes( 1, 1 ) = 2.0;
    green_edge_nodes( 1, 2 ) = 0.0;

    green_edge_node_normals( 0, 0 ) = 0.0;
    green_edge_node_normals( 0, 1 ) = 0.0;
    green_edge_node_normals( 0, 2 ) = 1.0;
    green_edge_node_normals( 1, 0 ) = 0.0;
    green_edge_node_normals( 1, 1 ) = 0.0;
    green_edge_node_normals( 1, 2 ) = 1.0;

    has_projection = PP::projectPointFeatureToEdgeFeature(
        parameters, blue_point, blue_normal, green_edge_nodes,
        green_edge_node_normals, projected_blue_point, green_node_id );

    TEST_ASSERT( has_projection );
    TEST_EQUALITY( green_node_id, 0 );
    TEST_EQUALITY( projected_blue_point( 0 ), 0.0 );
    TEST_EQUALITY( projected_blue_point( 1 ), 0.0 );
    TEST_EQUALITY( projected_blue_point( 2 ), 0.0 );

    // Projection 6.
    blue_point( 0 ) = 0.0;
    blue_point( 1 ) = 0.0;
    blue_point( 2 ) = 1.0;

    blue_normal( 0 ) = 0.0;
    blue_normal( 1 ) = 0.0;
    blue_normal( 2 ) = -1.0;

    green_edge_nodes( 0, 0 ) = 0.0;
    green_edge_nodes( 0, 1 ) = 0.0;
    green_edge_nodes( 0, 2 ) = 1.0;
    green_edge_nodes( 1, 0 ) = 2.0;
    green_edge_nodes( 1, 1 ) = 2.0;
    green_edge_nodes( 1, 2 ) = 1.0;

    green_edge_node_normals( 0, 0 ) = 0.0;
    green_edge_node_normals( 0, 1 ) = 0.0;
    green_edge_node_normals( 0, 2 ) = 1.0;
    green_edge_node_normals( 1, 0 ) = 0.0;
    green_edge_node_normals( 1, 1 ) = 0.0;
    green_edge_node_normals( 1, 2 ) = 1.0;

    has_projection = PP::projectPointFeatureToEdgeFeature(
        parameters, blue_point, blue_normal, green_edge_nodes,
        green_edge_node_normals, projected_blue_point, green_node_id );

    TEST_ASSERT( has_projection );
    TEST_EQUALITY( green_node_id, 0 );
    TEST_EQUALITY( projected_blue_point( 0 ), 0.0 );
    TEST_EQUALITY( projected_blue_point( 1 ), 0.0 );
    TEST_EQUALITY( projected_blue_point( 2 ), 1.0 );

    // Projection 7.
    blue_point( 0 ) = 1.0 + epsilon;
    blue_point( 1 ) = 1.0 + epsilon;
    blue_point( 2 ) = 1.0;

    blue_normal( 0 ) = 0.0;
    blue_normal( 1 ) = 0.0;
    blue_normal( 2 ) = -1.0;

    green_edge_nodes( 0, 0 ) = 0.0;
    green_edge_nodes( 0, 1 ) = 0.0;
    green_edge_nodes( 0, 2 ) = 0.0;
    green_edge_nodes( 1, 0 ) = 2.0;
    green_edge_nodes( 1, 1 ) = 2.0;
    green_edge_nodes( 1, 2 ) = 0.0;

    green_edge_node_normals( 0, 0 ) = 0.0;
    green_edge_node_normals( 0, 1 ) = 0.0;
    green_edge_node_normals( 0, 2 ) = 1.0;
    green_edge_node_normals( 1, 0 ) = 0.0;
    green_edge_node_normals( 1, 1 ) = 0.0;
    green_edge_node_normals( 1, 2 ) = 1.0;

    has_projection = PP::projectPointFeatureToEdgeFeature(
        parameters, blue_point, blue_normal, green_edge_nodes,
        green_edge_node_normals, projected_blue_point, green_node_id );

    TEST_ASSERT( has_projection );
    TEST_EQUALITY( green_node_id, -1 );
    TEST_FLOATING_EQUALITY( projected_blue_point( 0 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( projected_blue_point( 1 ), 1.0, epsilon );
    TEST_EQUALITY( projected_blue_point( 2 ), 0.0 );

    // Projection 8.
    blue_point( 0 ) = 1.0 + epsilon;
    blue_point( 1 ) = 1.0 - epsilon;
    blue_point( 2 ) = 1.0;

    blue_normal( 0 ) = 0.0;
    blue_normal( 1 ) = 0.0;
    blue_normal( 2 ) = -1.0;

    green_edge_nodes( 0, 0 ) = 0.0;
    green_edge_nodes( 0, 1 ) = 0.0;
    green_edge_nodes( 0, 2 ) = 0.0;
    green_edge_nodes( 1, 0 ) = 2.0;
    green_edge_nodes( 1, 1 ) = 2.0;
    green_edge_nodes( 1, 2 ) = 0.0;

    green_edge_node_normals( 0, 0 ) = 0.0;
    green_edge_node_normals( 0, 1 ) = 0.0;
    green_edge_node_normals( 0, 2 ) = 1.0;
    green_edge_node_normals( 1, 0 ) = 0.0;
    green_edge_node_normals( 1, 1 ) = 0.0;
    green_edge_node_normals( 1, 2 ) = 1.0;

    has_projection = PP::projectPointFeatureToEdgeFeature(
        parameters, blue_point, blue_normal, green_edge_nodes,
        green_edge_node_normals, projected_blue_point, green_node_id );

    TEST_ASSERT( has_projection );
    TEST_EQUALITY( green_node_id, -1 );
    TEST_FLOATING_EQUALITY( projected_blue_point( 0 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( projected_blue_point( 1 ), 1.0, epsilon );
    TEST_EQUALITY( projected_blue_point( 2 ), 0.0 );

    // Projection 9.
    blue_point( 0 ) = 1.0 - epsilon;
    blue_point( 1 ) = 1.0 + epsilon;
    blue_point( 2 ) = 1.0;

    blue_normal( 0 ) = 0.0;
    blue_normal( 1 ) = 0.0;
    blue_normal( 2 ) = -1.0;

    green_edge_nodes( 0, 0 ) = 0.0;
    green_edge_nodes( 0, 1 ) = 0.0;
    green_edge_nodes( 0, 2 ) = 0.0;
    green_edge_nodes( 1, 0 ) = 2.0;
    green_edge_nodes( 1, 1 ) = 2.0;
    green_edge_nodes( 1, 2 ) = 0.0;

    green_edge_node_normals( 0, 0 ) = 0.0;
    green_edge_node_normals( 0, 1 ) = 0.0;
    green_edge_node_normals( 0, 2 ) = 1.0;
    green_edge_node_normals( 1, 0 ) = 0.0;
    green_edge_node_normals( 1, 1 ) = 0.0;
    green_edge_node_normals( 1, 2 ) = 1.0;

    has_projection = PP::projectPointFeatureToEdgeFeature(
        parameters, blue_point, blue_normal, green_edge_nodes,
        green_edge_node_normals, projected_blue_point, green_node_id );

    TEST_ASSERT( has_projection );
    TEST_EQUALITY( green_node_id, -1 );
    TEST_FLOATING_EQUALITY( projected_blue_point( 0 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( projected_blue_point( 1 ), 1.0, epsilon );
    TEST_EQUALITY( projected_blue_point( 2 ), 0.0 );

    // Projection 10.
    blue_point( 0 ) = -epsilon;
    blue_point( 1 ) = epsilon;
    blue_point( 2 ) = 1.0;

    blue_normal( 0 ) = 0.0;
    blue_normal( 1 ) = 0.0;
    blue_normal( 2 ) = -1.0;

    green_edge_nodes( 0, 0 ) = 0.0;
    green_edge_nodes( 0, 1 ) = 0.0;
    green_edge_nodes( 0, 2 ) = 0.0;
    green_edge_nodes( 1, 0 ) = 2.0;
    green_edge_nodes( 1, 1 ) = 2.0;
    green_edge_nodes( 1, 2 ) = 0.0;

    green_edge_node_normals( 0, 0 ) = 0.0;
    green_edge_node_normals( 0, 1 ) = 0.0;
    green_edge_node_normals( 0, 2 ) = 1.0;
    green_edge_node_normals( 1, 0 ) = 0.0;
    green_edge_node_normals( 1, 1 ) = 0.0;
    green_edge_node_normals( 1, 2 ) = 1.0;

    has_projection = PP::projectPointFeatureToEdgeFeature(
        parameters, blue_point, blue_normal, green_edge_nodes,
        green_edge_node_normals, projected_blue_point, green_node_id );

    TEST_ASSERT( has_projection );
    TEST_EQUALITY( green_node_id, 0 );
    TEST_FLOATING_EQUALITY( projected_blue_point( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( projected_blue_point( 1 ), 0.0, epsilon );
    TEST_EQUALITY( projected_blue_point( 2 ), 0.0 );

    // Projection 11.
    blue_point( 0 ) = 4.80000000000000012e-04;
    blue_point( 1 ) = 3.48600000000000021e-02;
    blue_point( 2 ) = 7.00000000000000067e-02;

    blue_normal( 0 ) = -9.89672189476683117e-01;
    blue_normal( 1 ) = 0.00000000000000000e+00;
    blue_normal( 2 ) = -1.43349075254876335e-01;

    green_edge_nodes( 0, 0 ) = 4.79999999999980442e-04;
    green_edge_nodes( 0, 1 ) = 3.73199999999999990e-02;
    green_edge_nodes( 0, 2 ) = 6.94520547945205458e-02;
    green_edge_nodes( 1, 0 ) = 4.79999999999980442e-04;
    green_edge_nodes( 1, 1 ) = 3.73199999999999990e-02;
    green_edge_nodes( 1, 2 ) = 7.00000000000000067e-02;

    green_edge_node_normals( 0, 0 ) = 8.14680846865835639e-01;
    green_edge_node_normals( 0, 1 ) = -5.79909577218694516e-01;
    green_edge_node_normals( 0, 2 ) = 0.00000000000000000e+00;
    green_edge_node_normals( 1, 0 ) = 8.14680846865835639e-01;
    green_edge_node_normals( 1, 1 ) = -5.79909577218694516e-01;
    green_edge_node_normals( 1, 2 ) = 0.00000000000000000e+00;

    has_projection = PP::projectPointFeatureToEdgeFeature(
        parameters, blue_point, blue_normal, green_edge_nodes,
        green_edge_node_normals, projected_blue_point, green_node_id );

    TEST_ASSERT( has_projection );
    TEST_EQUALITY( green_node_id, 1 );
    TEST_FLOATING_EQUALITY( projected_blue_point( 0 ), green_edge_nodes( 1, 0 ),
                            epsilon );
    TEST_FLOATING_EQUALITY( projected_blue_point( 1 ), green_edge_nodes( 1, 1 ),
                            epsilon );
    TEST_FLOATING_EQUALITY( projected_blue_point( 2 ), green_edge_nodes( 1, 2 ),
                            epsilon );
}

//---------------------------------------------------------------------------//
// These tests have identical edge node normals creating a linear intersection
// problem.
TEUCHOS_UNIT_TEST( ProjectionPrimitives, linear_edge_edge_intersection_test )
{
    typedef DataTransferKit::ProjectionPrimitives PP;

    Teuchos::ParameterList parameters = defaultParameters();
    parameters.set<double>( "Merge Tolerance", epsilon );

    // Allocate Arrays.
    int space_dim = 3;
    Intrepid::FieldContainer<double> edge_1( 2, space_dim );
    Intrepid::FieldContainer<double> edge_2( 2, space_dim );
    Intrepid::FieldContainer<double> edge_2_node_normals( 2, space_dim );
    Intrepid::FieldContainer<double> edge_1_intersection( space_dim );
    Intrepid::FieldContainer<double> edge_2_intersection( space_dim );
    bool has_intersection = false;
    int node_id_1;
    int node_id_2;

    // Intersection 1. The edges are perpendicular in the xy plane and
    // separated in the z direction.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 1.0;
    edge_2( 0, 1 ) = 0.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 1.0;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, -1 );
    TEST_EQUALITY( node_id_2, -1 );

    TEST_FLOATING_EQUALITY( edge_1_intersection( 0 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 2 ), 0.0, epsilon );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 2 ), 1.0, epsilon );

    // Intersection 2. The edges are perpendicular in the xy plane and
    // in the same z plane.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 1.0;
    edge_2( 0, 1 ) = 0.0;
    edge_2( 0, 2 ) = 0.0;
    edge_2( 1, 0 ) = 1.0;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 0.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, -1 );
    TEST_EQUALITY( node_id_2, -1 );

    TEST_FLOATING_EQUALITY( edge_1_intersection( 0 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 2 ), 0.0, epsilon );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 2 ), 0.0, epsilon );

    // Intersection 3. The edges intersect at a the first blue node and the
    // first green node.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 0.0;
    edge_2( 0, 1 ) = 1.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 0.0;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 0 );
    TEST_EQUALITY( node_id_2, 0 );

    TEST_FLOATING_EQUALITY( edge_1_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 2 ), 0.0, epsilon );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 2 ), 1.0, epsilon );

    // Intersection 4. The edges intersect at the second blue node and the
    // first green node.
    edge_1( 0, 0 ) = 2.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 0.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 0.0;
    edge_2( 0, 1 ) = 1.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 0.0;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 1 );
    TEST_EQUALITY( node_id_2, 0 );

    TEST_FLOATING_EQUALITY( edge_1_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 2 ), 0.0, epsilon );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 2 ), 1.0, epsilon );

    // Intersection 5. The edges intersect at the first blue node and the
    // second green node.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 0.0;
    edge_2( 0, 1 ) = 2.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 0.0;
    edge_2( 1, 1 ) = 1.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 0 );
    TEST_EQUALITY( node_id_2, 1 );

    TEST_FLOATING_EQUALITY( edge_1_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 2 ), 0.0, epsilon );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 2 ), 1.0, epsilon );

    // Intersection 6. The edges intersect at the second blue node and the
    // second green node.
    edge_1( 0, 0 ) = 2.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 0.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 0.0;
    edge_2( 0, 1 ) = 2.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 0.0;
    edge_2( 1, 1 ) = 1.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 1 );
    TEST_EQUALITY( node_id_2, 1 );

    TEST_FLOATING_EQUALITY( edge_1_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 2 ), 0.0, epsilon );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 2 ), 1.0, epsilon );

    // Intersection 7. The edges intersect at the first blue node.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 0.0;
    edge_2( 0, 1 ) = 2.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 0.0;
    edge_2( 1, 1 ) = 0.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 0 );
    TEST_EQUALITY( node_id_2, -1 );

    TEST_FLOATING_EQUALITY( edge_1_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 2 ), 0.0, epsilon );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 2 ), 1.0, epsilon );

    // Intersection 8. The edges intersect at the first green node.
    edge_1( 0, 0 ) = -1.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 0.0;
    edge_2( 0, 1 ) = 1.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 0.0;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, -1 );
    TEST_EQUALITY( node_id_2, 0 );

    TEST_FLOATING_EQUALITY( edge_1_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 2 ), 0.0, epsilon );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 2 ), 1.0, epsilon );

    // Intersection 9. The edges intersect at the second blue node.
    edge_1( 0, 0 ) = 2.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 0.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 0.0;
    edge_2( 0, 1 ) = 0.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 0.0;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 1 );
    TEST_EQUALITY( node_id_2, -1 );

    TEST_FLOATING_EQUALITY( edge_1_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 2 ), 0.0, epsilon );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 2 ), 1.0, epsilon );

    // Intersection 10. The edges intersect at the second green node.
    edge_1( 0, 0 ) = 2.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = -1.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 0.0;
    edge_2( 0, 1 ) = 2.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 0.0;
    edge_2( 1, 1 ) = 1.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, -1 );
    TEST_EQUALITY( node_id_2, 1 );

    TEST_FLOATING_EQUALITY( edge_1_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 2 ), 0.0, epsilon );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 2 ), 1.0, epsilon );

    // Intersection 11. The edges do not overlap in the xy plane and therefore
    // do not intersect.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 3.0;
    edge_2( 0, 1 ) = 0.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 3.0;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );
    TEST_ASSERT( !has_intersection );

    // Intersection 12. The edges are the same, just shifted in z. We expect to
    // find the intersection at the first edge node.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 0.0;
    edge_2( 0, 1 ) = 1.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 0.0;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 0 );
    TEST_EQUALITY( node_id_2, 0 );

    TEST_FLOATING_EQUALITY( edge_1_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 2 ), 0.0, epsilon );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 0.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 1 ), 1.0, epsilon );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 2 ), 1.0, epsilon );

    // Intersection 13. The edges do not overlap in the xy plane and therefore
    // do not intersect. This time they differ by a small value close to the
    // floating point tolerance.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    double epsilon_scale = 100.0;
    edge_2( 0, 0 ) = 2.0 + epsilon_scale * epsilon;
    edge_2( 0, 1 ) = 0.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 2.0 + epsilon_scale * epsilon;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( !has_intersection );

    // Intersection 14. The edges overlap in the xy plane and therefore
    // intersect. This time they differ by a small value close to the
    // floating point tolerance.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 2.0 - epsilon_scale * epsilon;
    edge_2( 0, 1 ) = 0.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 2.0 - epsilon_scale * epsilon;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, -1 );
    TEST_EQUALITY( node_id_2, -1 );

    TEST_FLOATING_EQUALITY( edge_1_intersection( 0 ), 2.0,
                            epsilon_scale * epsilon );
    TEST_EQUALITY( edge_1_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_1_intersection( 2 ), 0.0 );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 2.0,
                            epsilon_scale * epsilon );
    TEST_EQUALITY( edge_2_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 2 ), 1.0 );

    // Intersection 15. The edges are identical, just shifted in z. This
    // should trigger a no intersection because because of ill-conditioning.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 0.0;
    edge_2( 0, 1 ) = 1.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 2.0;
    edge_2( 1, 1 ) = 1.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );
    TEST_ASSERT( !has_intersection );

    // Intersection 16. The edges overlap in the xy plane and therefore
    // intersect. This time they differ by a small value close to the
    // floating point tolerance such that the intersection point will be
    // shifted back to the original point. The epsilon perturbation is into
    // the x domain of the first edge.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 2.0 - 0.1 * epsilon;
    edge_2( 0, 1 ) = 0.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 2.0 - 0.1 * epsilon;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 1 );
    TEST_EQUALITY( node_id_2, -1 );

    TEST_EQUALITY( edge_1_intersection( 0 ), 2.0 );
    TEST_EQUALITY( edge_1_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_1_intersection( 2 ), 0.0 );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 2.0, epsilon );
    TEST_EQUALITY( edge_2_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 2 ), 1.0 );

    // Intersection 17. The edges overlap in the xy plane and therefore
    // intersect. This time they differ by a small value close to the
    // floating point tolerance such that the intersection point will be
    // shifted back to the original point. The epsilon perturbation is out of
    // the x domain of the first edge.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 2.0 + 0.1 * epsilon;
    edge_2( 0, 1 ) = 0.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 2.0 + 0.1 * epsilon;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 1 );
    TEST_EQUALITY( node_id_2, -1 );

    TEST_EQUALITY( edge_1_intersection( 0 ), 2.0 );
    TEST_EQUALITY( edge_1_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_1_intersection( 2 ), 0.0 );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 2.0, epsilon );
    TEST_EQUALITY( edge_2_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 2 ), 1.0 );

    // Intersection 16. The edges are identical with epsilon
    // perturbations. This should trigger a no intersection because because of
    // ill-conditioning.
    edge_1( 0, 0 ) = epsilon;
    edge_1( 0, 1 ) = 1.0 + epsilon;
    edge_1( 0, 2 ) = epsilon;
    edge_1( 1, 0 ) = 2.0 + epsilon;
    edge_1( 1, 1 ) = 1.0 + epsilon;
    edge_1( 1, 2 ) = epsilon;

    edge_2( 0, 0 ) = epsilon;
    edge_2( 0, 1 ) = 1.0 - epsilon;
    edge_2( 0, 2 ) = epsilon;
    edge_2( 1, 0 ) = 2.0;
    edge_2( 1, 1 ) = 1.0;
    edge_2( 1, 2 ) = 1.0 - epsilon;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = 0.0;
    edge_2_node_normals( 0, 2 ) = -1.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = 0.0;
    edge_2_node_normals( 1, 2 ) = -1.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );
    TEST_ASSERT( !has_intersection );
}

//---------------------------------------------------------------------------//
// These tests have different edge node normals that create a nonlinear
// intersection problem.
TEUCHOS_UNIT_TEST( ProjectionPrimitives, nonlinear_edge_edge_intersection_test )
{
    typedef DataTransferKit::ProjectionPrimitives PP;

    Teuchos::ParameterList parameters = defaultParameters();
    parameters.set<double>( "Merge Tolerance", epsilon );

    // Allocate Arrays.
    int space_dim = 3;
    Intrepid::FieldContainer<double> edge_1( 2, space_dim );
    Intrepid::FieldContainer<double> edge_2( 2, space_dim );
    Intrepid::FieldContainer<double> edge_2_node_normals( 2, space_dim );
    Intrepid::FieldContainer<double> edge_1_intersection( space_dim );
    Intrepid::FieldContainer<double> edge_2_intersection( space_dim );
    int node_id_1;
    int node_id_2;
    bool has_intersection = false;

    // Intersection 1. The edges are perpendicular in the xy plane and
    // separated in the z direction.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 1.0;
    edge_2( 0, 1 ) = 0.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 1.0;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 0, 2 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 2 ) = -std::sqrt( 2.0 ) / 2.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, -1 );
    TEST_EQUALITY( node_id_2, -1 );

    TEST_EQUALITY( edge_1_intersection( 0 ), 1.0 );
    TEST_EQUALITY( edge_1_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_1_intersection( 2 ), 0.0 );

    TEST_EQUALITY( edge_2_intersection( 0 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 2 ), 1.0 );

    // Intersection 2. The edges are perpendicular in the xy plane and
    // in the same z plane.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 1.0;
    edge_2( 0, 1 ) = 0.0;
    edge_2( 0, 2 ) = 0.0;
    edge_2( 1, 0 ) = 1.0;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 0.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 0, 2 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 2 ) = -std::sqrt( 2.0 ) / 2.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, -1 );
    TEST_EQUALITY( node_id_2, -1 );

    TEST_EQUALITY( edge_1_intersection( 0 ), 1.0 );
    TEST_EQUALITY( edge_1_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_1_intersection( 2 ), 0.0 );

    TEST_EQUALITY( edge_2_intersection( 0 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 2 ), 0.0 );

    // Intersection 3. The edges intersect at the first blue edge node.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 0.0;
    edge_2( 0, 1 ) = 1.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 0.0;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 0, 2 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 2 ) = -std::sqrt( 2.0 ) / 2.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 0 );
    TEST_EQUALITY( node_id_2, -1 );

    TEST_EQUALITY( edge_1_intersection( 0 ), 0.0 );
    TEST_EQUALITY( edge_1_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_1_intersection( 2 ), 0.0 );

    TEST_EQUALITY( edge_2_intersection( 0 ), 0.0 );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 1 ), 4.0 / 3.0, epsilon );
    TEST_EQUALITY( edge_2_intersection( 2 ), 1.0 );

    // Intersection 4. We expect to find the intersection at the second blue
    // edge node.
    edge_1( 0, 0 ) = 2.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 0.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 0.0;
    edge_2( 0, 1 ) = 1.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 0.0;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 0, 2 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 2 ) = -std::sqrt( 2.0 ) / 2.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 1 );
    TEST_EQUALITY( node_id_2, -1 );

    TEST_EQUALITY( edge_1_intersection( 0 ), 0.0 );
    TEST_EQUALITY( edge_1_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_1_intersection( 2 ), 0.0 );

    TEST_EQUALITY( edge_2_intersection( 0 ), 0.0 );
    TEST_FLOATING_EQUALITY( edge_2_intersection( 1 ), 4.0 / 3.0, epsilon );
    TEST_EQUALITY( edge_2_intersection( 2 ), 1.0 );

    // Intersection 5. We expect to find the intersection at the first green
    // edge node.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 1.0;
    edge_1( 1, 0 ) = 0.0;
    edge_1( 1, 1 ) = 2.0;
    edge_1( 1, 2 ) = 1.0;

    edge_2( 0, 0 ) = 1.0;
    edge_2( 0, 1 ) = 0.5;
    edge_2( 0, 2 ) = 0.0;
    edge_2( 1, 0 ) = 2.0;
    edge_2( 1, 1 ) = 0.5;
    edge_2( 1, 2 ) = 0.0;

    edge_2_node_normals( 0, 0 ) = -std::sqrt( 3.0 ) / 3.0;
    edge_2_node_normals( 0, 1 ) = std::sqrt( 3.0 ) / 3.0;
    edge_2_node_normals( 0, 2 ) = std::sqrt( 3.0 ) / 3.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 2 ) = std::sqrt( 2.0 ) / 2.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, -1 );
    TEST_EQUALITY( node_id_2, 0 );

    TEST_EQUALITY( edge_1_intersection( 0 ), 0.0 );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 1 ), 1.5, epsilon );
    TEST_EQUALITY( edge_1_intersection( 2 ), 1.0 );

    TEST_EQUALITY( edge_2_intersection( 0 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 1 ), 0.5 );
    TEST_EQUALITY( edge_2_intersection( 2 ), 0.0 );

    // Intersection 6. We expect to find the intersection at the second green
    // edge node.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 1.0;
    edge_1( 1, 0 ) = 0.0;
    edge_1( 1, 1 ) = 2.0;
    edge_1( 1, 2 ) = 1.0;

    edge_2( 0, 0 ) = 2.0;
    edge_2( 0, 1 ) = 0.5;
    edge_2( 0, 2 ) = 0.0;
    edge_2( 1, 0 ) = 1.0;
    edge_2( 1, 1 ) = 0.5;
    edge_2( 1, 2 ) = 0.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 0, 2 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 0 ) = -std::sqrt( 3.0 ) / 3.0;
    edge_2_node_normals( 1, 1 ) = std::sqrt( 3.0 ) / 3.0;
    edge_2_node_normals( 1, 2 ) = std::sqrt( 3.0 ) / 3.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, -1 );
    TEST_EQUALITY( node_id_2, 1 );

    TEST_EQUALITY( edge_1_intersection( 0 ), 0.0 );
    TEST_FLOATING_EQUALITY( edge_1_intersection( 1 ), 1.5, epsilon );
    TEST_EQUALITY( edge_1_intersection( 2 ), 1.0 );

    TEST_EQUALITY( edge_2_intersection( 0 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 1 ), 0.5 );
    TEST_EQUALITY( edge_2_intersection( 2 ), 0.0 );

    // Intersection 7. The edges intersect at the first blue edge node and
    // the first green edge node.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 2.0;
    edge_1( 0, 2 ) = 1.0;
    edge_1( 1, 0 ) = 0.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 1.0;

    edge_2( 0, 0 ) = 1.0;
    edge_2( 0, 1 ) = 1.0;
    edge_2( 0, 2 ) = 0.0;
    edge_2( 1, 0 ) = 2.0;
    edge_2( 1, 1 ) = 1.0;
    edge_2( 1, 2 ) = 0.0;

    edge_2_node_normals( 0, 0 ) = -std::sqrt( 3.0 ) / 3.0;
    edge_2_node_normals( 0, 1 ) = std::sqrt( 3.0 ) / 3.0;
    edge_2_node_normals( 0, 2 ) = std::sqrt( 3.0 ) / 3.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 2 ) = std::sqrt( 2.0 ) / 2.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 0 );
    TEST_EQUALITY( node_id_2, 0 );

    TEST_EQUALITY( edge_1_intersection( 0 ), 0.0 );
    TEST_EQUALITY( edge_1_intersection( 1 ), 2.0 );
    TEST_EQUALITY( edge_1_intersection( 2 ), 1.0 );

    TEST_EQUALITY( edge_2_intersection( 0 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 2 ), 0.0 );

    // Intersection 8. The edges intersect at the first blue edge node and
    // the second green edge node.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 2.0;
    edge_1( 0, 2 ) = 1.0;
    edge_1( 1, 0 ) = 0.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 1.0;

    edge_2( 0, 0 ) = 2.0;
    edge_2( 0, 1 ) = 1.0;
    edge_2( 0, 2 ) = 0.0;
    edge_2( 1, 0 ) = 1.0;
    edge_2( 1, 1 ) = 1.0;
    edge_2( 1, 2 ) = 0.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 0, 2 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 0 ) = -std::sqrt( 3.0 ) / 3.0;
    edge_2_node_normals( 1, 1 ) = std::sqrt( 3.0 ) / 3.0;
    edge_2_node_normals( 1, 2 ) = std::sqrt( 3.0 ) / 3.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 0 );
    TEST_EQUALITY( node_id_2, 1 );

    TEST_EQUALITY( edge_1_intersection( 0 ), 0.0 );
    TEST_EQUALITY( edge_1_intersection( 1 ), 2.0 );
    TEST_EQUALITY( edge_1_intersection( 2 ), 1.0 );

    TEST_EQUALITY( edge_2_intersection( 0 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 2 ), 0.0 );

    // Intersection 9. The edges intersect at the second blue edge node and
    // the first green edge node.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 1.0;
    edge_1( 1, 0 ) = 0.0;
    edge_1( 1, 1 ) = 2.0;
    edge_1( 1, 2 ) = 1.0;

    edge_2( 0, 0 ) = 1.0;
    edge_2( 0, 1 ) = 1.0;
    edge_2( 0, 2 ) = 0.0;
    edge_2( 1, 0 ) = 2.0;
    edge_2( 1, 1 ) = 1.0;
    edge_2( 1, 2 ) = 0.0;

    edge_2_node_normals( 0, 0 ) = -std::sqrt( 3.0 ) / 3.0;
    edge_2_node_normals( 0, 1 ) = std::sqrt( 3.0 ) / 3.0;
    edge_2_node_normals( 0, 2 ) = std::sqrt( 3.0 ) / 3.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 2 ) = std::sqrt( 2.0 ) / 2.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 1 );
    TEST_EQUALITY( node_id_2, 0 );

    TEST_EQUALITY( edge_1_intersection( 0 ), 0.0 );
    TEST_EQUALITY( edge_1_intersection( 1 ), 2.0 );
    TEST_EQUALITY( edge_1_intersection( 2 ), 1.0 );

    TEST_EQUALITY( edge_2_intersection( 0 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 2 ), 0.0 );

    // Intersection 10. The edges are the same, just shifted in z. We expect to
    // find the intersection at the second blue and green edge nodes.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 1.0;
    edge_1( 1, 0 ) = 0.0;
    edge_1( 1, 1 ) = 2.0;
    edge_1( 1, 2 ) = 1.0;

    edge_2( 0, 0 ) = 2.0;
    edge_2( 0, 1 ) = 1.0;
    edge_2( 0, 2 ) = 0.0;
    edge_2( 1, 0 ) = 0.0;
    edge_2( 1, 1 ) = 1.0;
    edge_2( 1, 2 ) = 0.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 0, 2 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 2 ) = std::sqrt( 2.0 ) / 2.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, 1 );
    TEST_EQUALITY( node_id_2, 1 );

    TEST_EQUALITY( edge_1_intersection( 0 ), 0.0 );
    TEST_EQUALITY( edge_1_intersection( 1 ), 2.0 );
    TEST_EQUALITY( edge_1_intersection( 2 ), 1.0 );

    TEST_EQUALITY( edge_2_intersection( 0 ), 0.0 );
    TEST_EQUALITY( edge_2_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 2 ), 0.0 );

    // Intersection 11. The edges do not overlap in the xy plane and therefore
    // do not intersect.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 3.0;
    edge_2( 0, 1 ) = 0.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 3.0;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 0, 2 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 2 ) = -std::sqrt( 2.0 ) / 2.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( !has_intersection );

    // Intersection 12. The edges do not overlap in the xy plane and therefore
    // do not intersect. This time they differ by a small value close to the
    // floating point tolerance.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    double epsilon_scale = 100.0;
    edge_2( 0, 0 ) = 2.0 + epsilon_scale * epsilon;
    edge_2( 0, 1 ) = 0.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 2.0 + epsilon_scale * epsilon;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 0, 2 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 2 ) = -std::sqrt( 2.0 ) / 2.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( !has_intersection );

    // Intersection 13. The edges overlap in the xy plane and therefore
    // intersect. This time they differ by a small value close to the
    // floating point tolerance.
    edge_1( 0, 0 ) = 0.0;
    edge_1( 0, 1 ) = 1.0;
    edge_1( 0, 2 ) = 0.0;
    edge_1( 1, 0 ) = 2.0;
    edge_1( 1, 1 ) = 1.0;
    edge_1( 1, 2 ) = 0.0;

    edge_2( 0, 0 ) = 2.0 - epsilon_scale * epsilon;
    edge_2( 0, 1 ) = 0.0;
    edge_2( 0, 2 ) = 1.0;
    edge_2( 1, 0 ) = 2.0 - epsilon_scale * epsilon;
    edge_2( 1, 1 ) = 2.0;
    edge_2( 1, 2 ) = 1.0;

    edge_2_node_normals( 0, 0 ) = 0.0;
    edge_2_node_normals( 0, 1 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 0, 2 ) = -std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 0 ) = 0.0;
    edge_2_node_normals( 1, 1 ) = std::sqrt( 2.0 ) / 2.0;
    edge_2_node_normals( 1, 2 ) = -std::sqrt( 2.0 ) / 2.0;

    has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( has_intersection );
    TEST_EQUALITY( node_id_1, -1 );
    TEST_EQUALITY( node_id_2, -1 );

    TEST_FLOATING_EQUALITY( edge_1_intersection( 0 ), 2.0,
                            epsilon_scale * epsilon );
    TEST_EQUALITY( edge_1_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_1_intersection( 2 ), 0.0 );

    TEST_FLOATING_EQUALITY( edge_2_intersection( 0 ), 2.0,
                            epsilon_scale * epsilon );
    TEST_EQUALITY( edge_2_intersection( 1 ), 1.0 );
    TEST_EQUALITY( edge_2_intersection( 2 ), 1.0 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ProjectionPrimitives, edge_edge_intersection_test_2 )
{
    typedef DataTransferKit::ProjectionPrimitives PP;

    Teuchos::ParameterList parameters = defaultParameters();
    parameters.set<double>( "Merge Tolerance", epsilon );

    int space_dim = 3;
    Intrepid::FieldContainer<double> edge_1( 2, space_dim );
    Intrepid::FieldContainer<double> edge_2( 2, space_dim );
    Intrepid::FieldContainer<double> edge_2_node_normals( 2, space_dim );
    Intrepid::FieldContainer<double> edge_1_intersection( space_dim );
    Intrepid::FieldContainer<double> edge_2_intersection( space_dim );
    int node_id_1;
    int node_id_2;

    edge_1( 0, 0 ) = 2.991988424774516;
    edge_1( 0, 1 ) = 0.0009914748115963934;
    edge_1( 0, 2 ) = 3.182025927246416;
    edge_1( 1, 0 ) = 2.967311343226939;
    edge_1( 1, 1 ) = 0.0005911593836794732;
    edge_1( 1, 2 ) = 4.092390551358848;

    edge_2( 0, 0 ) = 2.853291437028165;
    edge_2( 0, 1 ) = 0.9266757660485945;
    edge_2( 0, 2 ) = 3.361587438549589;
    edge_2( 1, 0 ) = 2.999999988912588;
    edge_2( 1, 1 ) = -0.0002579233771464424;
    edge_2( 1, 2 ) = 3.362794808389883;

    edge_2_node_normals( 0, 0 ) = 0.8902012057321598;
    edge_2_node_normals( 0, 1 ) = 0.4533729205081594;
    edge_2_node_normals( 0, 2 ) = -0.04466327644621499;
    edge_2_node_normals( 1, 0 ) = 0.9866969130275373;
    edge_2_node_normals( 1, 1 ) = -0.1562983090786018;
    edge_2_node_normals( 1, 2 ) = -0.04472181124572589;

    bool has_intersection = PP::edgeEdgeIntersection(
        parameters, edge_1, edge_2, edge_2_node_normals, edge_1_intersection,
        edge_2_intersection, node_id_1, node_id_2 );

    TEST_ASSERT( !has_intersection );
}

//---------------------------------------------------------------------------//
// end tstProjectionPrimitives.cpp
//---------------------------------------------------------------------------//
