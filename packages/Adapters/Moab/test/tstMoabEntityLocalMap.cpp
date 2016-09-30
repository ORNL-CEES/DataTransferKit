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
 * \file tstMoabEntityLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief MoabEntityLocalMap unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_MoabEntity.hpp>
#include <DTK_MoabEntityExtraData.hpp>
#include <DTK_MoabMeshSetIndexer.hpp>
#include <DTK_MoabEntityLocalMap.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

#include <moab/Interface.hpp>
#include <moab/ParallelComm.hpp>
#include <moab/Core.hpp>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// TEST EPSILON
//---------------------------------------------------------------------------//

const double epsilon = 1.0e-14;

//---------------------------------------------------------------------------//
// Hex-8 test.
TEUCHOS_UNIT_TEST( MoabEntity, hex_8_test )
{
    // Extract the raw mpi communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();
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
    Teuchos::Array<moab::EntityHandle> nodes(8);
    double node_coords[3];
    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = -2.0;
    error = moab_mesh->create_vertex( node_coords, nodes[0] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 2.0;
    node_coords[1] = 0.0;
    node_coords[2] = -2.0;
    error = moab_mesh->create_vertex( node_coords, nodes[1] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 2.0;
    node_coords[1] = 2.0;
    node_coords[2] = -2.0;
    error = moab_mesh->create_vertex( node_coords, nodes[2] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 0.0;
    node_coords[1] = 2.0;
    node_coords[2] = -2.0;
    error = moab_mesh->create_vertex( node_coords, nodes[3] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 0.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;
    error = moab_mesh->create_vertex( node_coords, nodes[4] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 2.0;
    node_coords[1] = 0.0;
    node_coords[2] = 0.0;
    error = moab_mesh->create_vertex( node_coords, nodes[5] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 2.0;
    node_coords[1] = 2.0;
    node_coords[2] = 0.0;
    error = moab_mesh->create_vertex( node_coords, nodes[6] );
    TEST_EQUALITY( error, moab::MB_SUCCESS );

    node_coords[0] = 0.0;
    node_coords[1] = 2.0;
    node_coords[2] = 0.0;
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

    // Create a local map from the moab mesh.
    Teuchos::RCP<DataTransferKit::EntityLocalMap> local_map =
        Teuchos::rcp( new DataTransferKit::MoabEntityLocalMap(parallel_mesh) );

    // Test the measure.
    TEST_EQUALITY( local_map->measure(dtk_entity), 8.0 );

    // Test the centroid.
    Teuchos::Array<double> centroid( space_dim, 0.0 );
    local_map->centroid( dtk_entity, centroid() );
    TEST_EQUALITY( centroid[0], 1.0 );
    TEST_EQUALITY( centroid[1], 1.0 );
    TEST_EQUALITY( centroid[2], -1.0 );

    // Make a good point and a bad point.
    Teuchos::Array<double> good_point( space_dim );
    good_point[0] = 0.5;
    good_point[1] = 1.5;
    good_point[2] = -1.0;
    Teuchos::Array<double> bad_point( space_dim );
    bad_point[0] = 0.75;
    bad_point[1] = -1.75;
    bad_point[2] = 0.35;

    // Test the reference frame safeguard.
    TEST_ASSERT(
            local_map->isSafeToMapToReferenceFrame(dtk_entity,good_point()) );
    TEST_ASSERT(
            !local_map->isSafeToMapToReferenceFrame(dtk_entity,bad_point()) );

    // Test the mapping to reference frame.
    Teuchos::Array<double> ref_good_point( space_dim );
    bool good_map = local_map->mapToReferenceFrame(
            dtk_entity, good_point(), ref_good_point() );
    TEST_ASSERT( good_map );
    TEST_FLOATING_EQUALITY( ref_good_point[0], -0.5, epsilon );
    TEST_FLOATING_EQUALITY( ref_good_point[1], 0.5, epsilon );
    TEST_ASSERT( std::abs(ref_good_point[2]) < epsilon );

    Teuchos::Array<double> ref_bad_point( space_dim );
    bool bad_map = local_map->mapToReferenceFrame(
            dtk_entity, bad_point(), ref_bad_point() );
    TEST_ASSERT( !bad_map );

    // Test the point inclusion.
    TEST_ASSERT( local_map->checkPointInclusion(dtk_entity,ref_good_point()) );
    TEST_ASSERT( !local_map->checkPointInclusion(dtk_entity,ref_bad_point()) );

    // Test the map to physical frame.
    Teuchos::Array<double> phy_good_point( space_dim );
    local_map->mapToPhysicalFrame(dtk_entity,ref_good_point(),phy_good_point());
    TEST_FLOATING_EQUALITY( good_point[0], phy_good_point[0], epsilon );
    TEST_FLOATING_EQUALITY( good_point[1], phy_good_point[1], epsilon );
    TEST_FLOATING_EQUALITY( good_point[2], phy_good_point[2], epsilon );

    Teuchos::Array<double> phy_bad_point( space_dim );
    local_map->mapToPhysicalFrame(dtk_entity,ref_bad_point(),phy_bad_point());
    TEST_FLOATING_EQUALITY( bad_point[0], phy_bad_point[0], epsilon );
    TEST_FLOATING_EQUALITY( bad_point[1], phy_bad_point[1], epsilon );
    TEST_FLOATING_EQUALITY( bad_point[2], phy_bad_point[2], epsilon );

    // Test the coordinates of the points extracted through the centroid
    // function.
    DataTransferKit::Entity dtk_node;
    Teuchos::Array<double> point_coords(space_dim);
    int num_nodes = 8;
    for ( int n = 0; n < num_nodes; ++n )
    {
        dtk_node = DataTransferKit::MoabEntity(
            nodes[n], parallel_mesh.ptr(), set_indexer.ptr() );
        local_map->centroid( dtk_node, point_coords() );
        moab_mesh->get_coords( &nodes[n], 1, node_coords );
        TEST_EQUALITY( node_coords[0], point_coords[0] );
        TEST_EQUALITY( node_coords[1], point_coords[1] );
        TEST_EQUALITY( node_coords[2], point_coords[2] );
    }
}

//---------------------------------------------------------------------------//
// Quad-4 test.
TEUCHOS_UNIT_TEST( MoabEntity, quad_4_test )
{
    // Extract the raw mpi communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();
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
    Teuchos::Array<moab::EntityHandle> nodes(4);
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

    // Create a local map from the moab mesh.
    Teuchos::RCP<DataTransferKit::EntityLocalMap> local_map =
        Teuchos::rcp( new DataTransferKit::MoabEntityLocalMap(parallel_mesh) );

    // Test the measure.
    TEST_EQUALITY( local_map->measure(dtk_entity), 4.0 );

    // Test the centroid.
    Teuchos::Array<double> centroid( space_dim, 0.0 );
    local_map->centroid( dtk_entity, centroid() );
    TEST_EQUALITY( centroid[0], 1.0 );
    TEST_EQUALITY( centroid[1], 1.0 );

    // Make a good point and a bad point.
    Teuchos::Array<double> good_point( space_dim );
    good_point[0] = 0.5;
    good_point[1] = 1.5;
    Teuchos::Array<double> bad_point( space_dim );
    bad_point[0] = 0.75;
    bad_point[1] = -1.75;

    // Test the reference frame safeguard.
    TEST_ASSERT(
            local_map->isSafeToMapToReferenceFrame(dtk_entity,good_point()) );
    TEST_ASSERT(
            !local_map->isSafeToMapToReferenceFrame(dtk_entity,bad_point()) );

    // Test the mapping to reference frame.
    Teuchos::Array<double> ref_good_point( space_dim );
    bool good_map = local_map->mapToReferenceFrame(
            dtk_entity, good_point(), ref_good_point() );
    TEST_ASSERT( good_map );
    TEST_FLOATING_EQUALITY( ref_good_point[0], -0.5, epsilon );
    TEST_FLOATING_EQUALITY( ref_good_point[1], 0.5, epsilon );

    Teuchos::Array<double> ref_bad_point( space_dim );
    bool bad_map = local_map->mapToReferenceFrame(
            dtk_entity, bad_point(), ref_bad_point() );
    TEST_ASSERT( !bad_map );

    // Test the point inclusion.
    TEST_ASSERT( local_map->checkPointInclusion(dtk_entity,ref_good_point()) );
    TEST_ASSERT( !local_map->checkPointInclusion(dtk_entity,ref_bad_point()) );

    // Test the map to physical frame.
    Teuchos::Array<double> phy_good_point( space_dim );
    local_map->mapToPhysicalFrame(dtk_entity,ref_good_point(),phy_good_point());
    TEST_FLOATING_EQUALITY( good_point[0], phy_good_point[0], epsilon );
    TEST_FLOATING_EQUALITY( good_point[1], phy_good_point[1], epsilon );

    Teuchos::Array<double> phy_bad_point( space_dim );
    local_map->mapToPhysicalFrame(dtk_entity,ref_bad_point(),phy_bad_point());
    TEST_FLOATING_EQUALITY( bad_point[0], phy_bad_point[0], epsilon );
    TEST_FLOATING_EQUALITY( bad_point[1], phy_bad_point[1], epsilon );

    // Test the coordinates of the points extracted through the centroid
    // function.
    DataTransferKit::Entity dtk_node;
    Teuchos::Array<double> point_coords(space_dim);
    int num_nodes = 4;
    for ( int n = 0; n < num_nodes; ++n )
    {
        dtk_node = DataTransferKit::MoabEntity(
            nodes[n], parallel_mesh.ptr(), set_indexer.ptr() );
        local_map->centroid( dtk_node, point_coords() );
        moab_mesh->get_coords( &nodes[n], 1, node_coords );
        TEST_EQUALITY( node_coords[0], point_coords[0] );
        TEST_EQUALITY( node_coords[1], point_coords[1] );
    }
}

//---------------------------------------------------------------------------//
// end tstMoabEntityLocalMap.cpp
//---------------------------------------------------------------------------//

