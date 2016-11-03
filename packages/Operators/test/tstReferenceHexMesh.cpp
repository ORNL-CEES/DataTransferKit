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
 * \file tstReferenceNode.cpp
 * \author Stuart R. Slattery
 * \brief ReferenceNode unit tests.
 */
//---------------------------------------------------------------------------//

#include "DTK_BasicEntitySet.hpp"
#include "DTK_Jacobian.hpp"
#include "reference_implementation/DTK_ReferenceHex.hpp"
#include "reference_implementation/DTK_ReferenceHexIntegrationRule.hpp"
#include "reference_implementation/DTK_ReferenceHexLocalMap.hpp"
#include "reference_implementation/DTK_ReferenceHexMesh.hpp"
#include "reference_implementation/DTK_ReferenceHexShapeFunction.hpp"
#include "reference_implementation/DTK_ReferenceNode.hpp"

#include <DTK_BasicEntityPredicates.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <limits>

//---------------------------------------------------------------------------//
// TEST EPSILON
//---------------------------------------------------------------------------//

const double epsilon = 100.0 * std::numeric_limits<double>::epsilon();

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ReferenceHexMesh, cell_constructor_test )
{
    // Get the comm.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Create the mesh.
    unsigned x_num_cells = 20;
    unsigned y_num_cells = 12;
    unsigned local_z_num_cells = 15;
    unsigned z_num_cells = local_z_num_cells * comm_size;

    unsigned local_num_cells = x_num_cells * y_num_cells * local_z_num_cells;
    unsigned global_num_cells = x_num_cells * y_num_cells * z_num_cells;

    unsigned x_num_nodes = x_num_cells + 1;
    unsigned y_num_nodes = y_num_cells + 1;
    unsigned local_z_num_nodes = local_z_num_cells + 1;
    unsigned z_num_nodes = z_num_cells + 1;

    unsigned local_num_nodes = x_num_nodes * y_num_nodes;
    if ( comm_rank < comm_size - 1 )
    {
        local_num_nodes *= local_z_num_nodes - 1;
    }
    else
    {
        local_num_nodes *= local_z_num_nodes;
    }

    unsigned global_num_nodes = x_num_nodes * y_num_nodes * z_num_nodes;

    double cell_width = 0.1;
    double x_min = 0.0;
    double x_max = cell_width * x_num_cells;
    double y_min = 0.0;
    double y_max = cell_width * y_num_cells;
    double z_min = 0.0;
    double z_max = cell_width * z_num_cells;

    DataTransferKit::UnitTest::ReferenceHexMesh mesh(
        comm, x_min, x_max, x_num_cells, y_min, y_max, y_num_cells, z_min,
        z_max, z_num_cells );

    // Test the entity set.
    auto entity_set = mesh.functionSpace()->entitySet();
    TEST_EQUALITY( entity_set->physicalDimension(), 3 );

    Teuchos::Tuple<double, 6> box;
    entity_set->localBoundingBox( box );
    TEST_EQUALITY( box[0], 0.0 );
    TEST_EQUALITY( box[1], 0.0 );
    TEST_EQUALITY( box[2], cell_width * local_z_num_cells * comm_rank );
    TEST_EQUALITY( box[3], cell_width * x_num_cells );
    TEST_EQUALITY( box[4], cell_width * y_num_cells );
    TEST_EQUALITY( box[5], cell_width * local_z_num_cells * ( comm_rank + 1 ) );

    entity_set->globalBoundingBox( box );
    TEST_EQUALITY( box[0], 0.0 );
    TEST_EQUALITY( box[1], 0.0 );
    TEST_EQUALITY( box[2], 0.0 );
    TEST_EQUALITY( box[3], cell_width * x_num_cells );
    TEST_EQUALITY( box[4], cell_width * y_num_cells );
    TEST_EQUALITY( box[5], cell_width * z_num_cells );

    // Check the nodes.
    DataTransferKit::Entity entity;
    for ( unsigned k = local_z_num_cells * comm_rank;
          k < local_z_num_cells * comm_rank + local_z_num_nodes; ++k )
    {
        for ( unsigned j = 0; j < y_num_nodes; ++j )
        {
            for ( unsigned i = 0; i < x_num_nodes; ++i )
            {
                DataTransferKit::EntityId node_id =
                    i + j * x_num_nodes + k * x_num_nodes * y_num_nodes;
                TEST_ASSERT( node_id < global_num_nodes );

                entity_set->getEntity( node_id, 0, entity );

                entity.boundingBox( box );
                TEST_FLOATING_EQUALITY( box[0], i * cell_width, epsilon );
                TEST_FLOATING_EQUALITY( box[1], j * cell_width, epsilon );
                TEST_FLOATING_EQUALITY( box[2], k * cell_width, epsilon );

                if ( k < local_z_num_cells * comm_rank + local_z_num_nodes -
                             1 ||
                     comm_rank == comm_size - 1 )
                {
                    TEST_EQUALITY( entity.ownerRank(), comm_rank );
                }
                else
                {
                    TEST_EQUALITY( entity.ownerRank(), comm_rank + 1 );
                }
            }
        }
    }

    // Check the cells.
    for ( unsigned k = local_z_num_cells * comm_rank;
          k < local_z_num_cells * ( comm_rank + 1 ); ++k )
    {
        for ( unsigned j = 0; j < y_num_cells; ++j )
        {
            for ( unsigned i = 0; i < x_num_cells; ++i )
            {
                DataTransferKit::EntityId cell_id =
                    i + j * x_num_cells + k * x_num_cells * y_num_cells;
                TEST_ASSERT( cell_id < global_num_cells );

                entity_set->getEntity( cell_id, 3, entity );

                entity.boundingBox( box );
                TEST_FLOATING_EQUALITY( box[0], i * cell_width, epsilon );
                TEST_FLOATING_EQUALITY( box[1], j * cell_width, epsilon );
                TEST_FLOATING_EQUALITY( box[2], k * cell_width, epsilon );
                TEST_FLOATING_EQUALITY( box[3], ( i + 1 ) * cell_width,
                                        epsilon );
                TEST_FLOATING_EQUALITY( box[4], ( j + 1 ) * cell_width,
                                        epsilon );
                TEST_FLOATING_EQUALITY( box[5], ( k + 1 ) * cell_width,
                                        epsilon );

                TEST_EQUALITY( entity.ownerRank(), comm_rank );
            }
        }
    }

    // Check iterators.
    DataTransferKit::LocalEntityPredicate pred( comm_rank );
    DataTransferKit::EntityIterator node_it =
        entity_set->entityIterator( 0, pred.getFunction() );
    TEST_EQUALITY( node_it.size(), local_num_nodes );
    DataTransferKit::EntityIterator cell_it =
        entity_set->entityIterator( 3, pred.getFunction() );
    TEST_EQUALITY( cell_it.size(), local_num_cells );

    // Check element construction via mapping and shape functions.
    auto local_map = mesh.functionSpace()->localMap();
    auto shape_function = mesh.functionSpace()->shapeFunction();

    // Nodes
    Teuchos::Array<DataTransferKit::SupportId> support_ids;
    Teuchos::Array<double> node_coords( 3 );
    auto nodes_begin = node_it.begin();
    auto nodes_end = node_it.end();
    for ( node_it = nodes_begin; node_it != nodes_end; ++node_it )
    {
        int k = node_it->id() / ( x_num_nodes * y_num_nodes );
        int j = ( node_it->id() - k * x_num_nodes * y_num_nodes ) / x_num_nodes;
        int i = node_it->id() - j * x_num_nodes - k * x_num_nodes * y_num_nodes;
        TEST_EQUALITY( node_it->id(),
                       i + j * x_num_nodes + k * x_num_nodes * y_num_nodes );

        TEST_EQUALITY( local_map->measure( *node_it ), 0.0 );

        local_map->centroid( *node_it, node_coords() );
        TEST_FLOATING_EQUALITY( node_coords[0], cell_width * i, epsilon );
        TEST_FLOATING_EQUALITY( node_coords[1], cell_width * j, epsilon );
        TEST_FLOATING_EQUALITY( node_coords[2], cell_width * k, epsilon );

        shape_function->entitySupportIds( *node_it, support_ids );
        TEST_EQUALITY( support_ids.size(), 1 );
        TEST_EQUALITY( support_ids[0], node_it->id() );
    }

    // Cells
    Teuchos::Array<double> centroid( 3 );
    Teuchos::Array<double> bad_point( 3, -0.1 );
    Teuchos::Array<double> ref_point( 3 );
    Teuchos::Array<double> bad_ref_point( 3, -1.1 );
    Teuchos::Array<double> phys_point( 3 );
    double volume = cell_width * cell_width * cell_width;
    auto cells_begin = cell_it.begin();
    auto cells_end = cell_it.end();
    for ( cell_it = cells_begin; cell_it != cells_end; ++cell_it )
    {
        int k = cell_it->id() / ( x_num_cells * y_num_cells );
        int j = ( cell_it->id() - k * x_num_cells * y_num_cells ) / x_num_cells;
        int i = cell_it->id() - j * x_num_cells - k * x_num_cells * y_num_cells;
        TEST_EQUALITY( cell_it->id(),
                       i + j * x_num_cells + k * x_num_cells * y_num_cells );

        TEST_FLOATING_EQUALITY( local_map->measure( *cell_it ), volume,
                                epsilon );

        local_map->centroid( *cell_it, centroid() );
        TEST_FLOATING_EQUALITY( centroid[0], cell_width * i + cell_width / 2,
                                epsilon );
        TEST_FLOATING_EQUALITY( centroid[1], cell_width * j + cell_width / 2,
                                epsilon );
        TEST_FLOATING_EQUALITY( centroid[2], cell_width * k + cell_width / 2,
                                epsilon );

        TEST_ASSERT(
            local_map->isSafeToMapToReferenceFrame( *cell_it, centroid() ) );
        TEST_ASSERT(
            !local_map->isSafeToMapToReferenceFrame( *cell_it, bad_point() ) );

        TEST_ASSERT( local_map->mapToReferenceFrame( *cell_it, centroid(),
                                                     ref_point() ) );
        TEST_EQUALITY( ref_point[0], 0.0 );
        TEST_EQUALITY( ref_point[1], 0.0 );
        TEST_EQUALITY( ref_point[2], 0.0 );

        TEST_ASSERT( local_map->checkPointInclusion( *cell_it, ref_point() ) );
        TEST_ASSERT(
            !local_map->checkPointInclusion( *cell_it, bad_ref_point() ) );

        local_map->mapToPhysicalFrame( *cell_it, ref_point(), phys_point() );
        TEST_EQUALITY( phys_point[0], centroid[0] );
        TEST_EQUALITY( phys_point[1], centroid[1] );
        TEST_EQUALITY( phys_point[2], centroid[2] );

        shape_function->entitySupportIds( *cell_it, support_ids );
        TEST_EQUALITY( support_ids.size(), 8 );
        TEST_EQUALITY( support_ids[0], ( i ) + (j)*x_num_nodes +
                                           (k)*x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[1], ( i + 1 ) + (j)*x_num_nodes +
                                           (k)*x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[2], ( i + 1 ) + ( j + 1 ) * x_num_nodes +
                                           (k)*x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[3], ( i ) + ( j + 1 ) * x_num_nodes +
                                           (k)*x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[4],
                       ( i ) + (j)*x_num_nodes +
                           ( k + 1 ) * x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[5],
                       ( i + 1 ) + (j)*x_num_nodes +
                           ( k + 1 ) * x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[6],
                       ( i + 1 ) + ( j + 1 ) * x_num_nodes +
                           ( k + 1 ) * x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[7],
                       ( i ) + ( j + 1 ) * x_num_nodes +
                           ( k + 1 ) * x_num_nodes * y_num_nodes );
    }

    // Check the nodal field.
    auto field = mesh.nodalField( 3 );
    auto field_supports = field->getLocalSupportIds();
    int counter = 0;
    for ( node_it = nodes_begin; node_it != nodes_end; ++node_it, ++counter )
    {
        TEST_EQUALITY( node_it->id(), field_supports[counter] );
    }
    for ( node_it = nodes_begin; node_it != nodes_end; ++node_it )
    {
        field->writeFieldData( node_it->id(), 0, node_it->id() + 32.2 );
        field->writeFieldData( node_it->id(), 1, node_it->id() + 3.2 );
    }

    for ( node_it = nodes_begin; node_it != nodes_end; ++node_it )
    {
        TEST_EQUALITY( field->readFieldData( node_it->id(), 0 ),
                       node_it->id() + 32.2 );
        TEST_EQUALITY( field->readFieldData( node_it->id(), 1 ),
                       node_it->id() + 3.2 );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ReferenceHexMesh, edge_constructor_test )
{
    // Get the comm.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Create the mesh.
    unsigned x_num_cells = 20;
    unsigned y_num_cells = 12;
    unsigned local_z_num_cells = 15;
    unsigned z_num_cells = local_z_num_cells * comm_size;

    unsigned local_num_cells = x_num_cells * y_num_cells * local_z_num_cells;
    unsigned global_num_cells = x_num_cells * y_num_cells * z_num_cells;

    unsigned x_num_nodes = x_num_cells + 1;
    unsigned y_num_nodes = y_num_cells + 1;
    unsigned local_z_num_nodes = local_z_num_cells + 1;
    unsigned z_num_nodes = z_num_cells + 1;

    unsigned local_num_nodes = x_num_nodes * y_num_nodes;
    if ( comm_rank < comm_size - 1 )
    {
        local_num_nodes *= local_z_num_nodes - 1;
    }
    else
    {
        local_num_nodes *= local_z_num_nodes;
    }

    unsigned global_num_nodes = x_num_nodes * y_num_nodes * z_num_nodes;

    Teuchos::Array<double> x_edges( x_num_nodes );
    Teuchos::Array<double> y_edges( y_num_nodes );
    Teuchos::Array<double> z_edges( z_num_nodes );

    double cell_width = 0.1;
    for ( unsigned i = 0; i < x_num_nodes; ++i )
        x_edges[i] = i * cell_width;
    for ( unsigned i = 0; i < y_num_nodes; ++i )
        y_edges[i] = i * cell_width;
    for ( unsigned i = 0; i < z_num_nodes; ++i )
        z_edges[i] = i * cell_width;

    DataTransferKit::UnitTest::ReferenceHexMesh mesh( comm, x_edges, y_edges,
                                                      z_edges );

    // Test the entity set.
    auto entity_set = mesh.functionSpace()->entitySet();
    TEST_EQUALITY( entity_set->physicalDimension(), 3 );

    Teuchos::Tuple<double, 6> box;
    entity_set->localBoundingBox( box );
    TEST_EQUALITY( box[0], 0.0 );
    TEST_EQUALITY( box[1], 0.0 );
    TEST_EQUALITY( box[2], cell_width * local_z_num_cells * comm_rank );
    TEST_EQUALITY( box[3], cell_width * x_num_cells );
    TEST_EQUALITY( box[4], cell_width * y_num_cells );
    TEST_EQUALITY( box[5], cell_width * local_z_num_cells * ( comm_rank + 1 ) );

    entity_set->globalBoundingBox( box );
    TEST_EQUALITY( box[0], 0.0 );
    TEST_EQUALITY( box[1], 0.0 );
    TEST_EQUALITY( box[2], 0.0 );
    TEST_EQUALITY( box[3], cell_width * x_num_cells );
    TEST_EQUALITY( box[4], cell_width * y_num_cells );
    TEST_EQUALITY( box[5], cell_width * z_num_cells );

    // Check the nodes.
    DataTransferKit::Entity entity;
    for ( unsigned k = local_z_num_cells * comm_rank;
          k < local_z_num_cells * comm_rank + local_z_num_nodes; ++k )
    {
        for ( unsigned j = 0; j < y_num_nodes; ++j )
        {
            for ( unsigned i = 0; i < x_num_nodes; ++i )
            {
                DataTransferKit::EntityId node_id =
                    i + j * x_num_nodes + k * x_num_nodes * y_num_nodes;
                TEST_ASSERT( node_id < global_num_nodes );

                entity_set->getEntity( node_id, 0, entity );

                entity.boundingBox( box );
                TEST_EQUALITY( box[0], i * cell_width );
                TEST_EQUALITY( box[1], j * cell_width );
                TEST_EQUALITY( box[2], k * cell_width );

                if ( k < local_z_num_cells * comm_rank + local_z_num_nodes -
                             1 ||
                     comm_rank == comm_size - 1 )
                {
                    TEST_EQUALITY( entity.ownerRank(), comm_rank );
                }
                else
                {
                    TEST_EQUALITY( entity.ownerRank(), comm_rank + 1 );
                }
            }
        }
    }

    // Check the cells.
    for ( unsigned k = local_z_num_cells * comm_rank;
          k < local_z_num_cells * ( comm_rank + 1 ); ++k )
    {
        for ( unsigned j = 0; j < y_num_cells; ++j )
        {
            for ( unsigned i = 0; i < x_num_cells; ++i )
            {
                DataTransferKit::EntityId cell_id =
                    i + j * x_num_cells + k * x_num_cells * y_num_cells;
                TEST_ASSERT( cell_id < global_num_cells );

                entity_set->getEntity( cell_id, 3, entity );

                entity.boundingBox( box );
                TEST_EQUALITY( box[0], i * cell_width );
                TEST_EQUALITY( box[1], j * cell_width );
                TEST_EQUALITY( box[2], k * cell_width );
                TEST_EQUALITY( box[3], ( i + 1 ) * cell_width );
                TEST_EQUALITY( box[4], ( j + 1 ) * cell_width );
                TEST_EQUALITY( box[5], ( k + 1 ) * cell_width );

                TEST_EQUALITY( entity.ownerRank(), comm_rank );
            }
        }
    }

    // Check iterators.
    DataTransferKit::LocalEntityPredicate pred( comm_rank );
    DataTransferKit::EntityIterator node_it =
        entity_set->entityIterator( 0, pred.getFunction() );
    TEST_EQUALITY( node_it.size(), local_num_nodes );
    DataTransferKit::EntityIterator cell_it =
        entity_set->entityIterator( 3, pred.getFunction() );
    TEST_EQUALITY( cell_it.size(), local_num_cells );

    // Check element construction via mapping and shape functions.
    auto local_map = mesh.functionSpace()->localMap();
    auto shape_function = mesh.functionSpace()->shapeFunction();

    // Nodes
    Teuchos::Array<DataTransferKit::SupportId> support_ids;
    Teuchos::Array<double> node_coords( 3 );
    auto nodes_begin = node_it.begin();
    auto nodes_end = node_it.end();
    for ( node_it = nodes_begin; node_it != nodes_end; ++node_it )
    {
        int k = node_it->id() / ( x_num_nodes * y_num_nodes );
        int j = ( node_it->id() - k * x_num_nodes * y_num_nodes ) / x_num_nodes;
        int i = node_it->id() - j * x_num_nodes - k * x_num_nodes * y_num_nodes;
        TEST_EQUALITY( node_it->id(),
                       i + j * x_num_nodes + k * x_num_nodes * y_num_nodes );

        TEST_EQUALITY( local_map->measure( *node_it ), 0.0 );

        local_map->centroid( *node_it, node_coords() );
        TEST_EQUALITY( node_coords[0], cell_width * i );
        TEST_EQUALITY( node_coords[1], cell_width * j );
        TEST_EQUALITY( node_coords[2], cell_width * k );

        shape_function->entitySupportIds( *node_it, support_ids );
        TEST_EQUALITY( support_ids.size(), 1 );
        TEST_EQUALITY( support_ids[0], node_it->id() );
    }

    // Cells
    Teuchos::Array<double> centroid( 3 );
    Teuchos::Array<double> bad_point( 3, -0.1 );
    Teuchos::Array<double> ref_point( 3 );
    Teuchos::Array<double> bad_ref_point( 3, -1.1 );
    Teuchos::Array<double> phys_point( 3 );
    double volume = cell_width * cell_width * cell_width;
    auto cells_begin = cell_it.begin();
    auto cells_end = cell_it.end();
    for ( cell_it = cells_begin; cell_it != cells_end; ++cell_it )
    {
        int k = cell_it->id() / ( x_num_cells * y_num_cells );
        int j = ( cell_it->id() - k * x_num_cells * y_num_cells ) / x_num_cells;
        int i = cell_it->id() - j * x_num_cells - k * x_num_cells * y_num_cells;
        TEST_EQUALITY( cell_it->id(),
                       i + j * x_num_cells + k * x_num_cells * y_num_cells );

        TEST_FLOATING_EQUALITY( local_map->measure( *cell_it ), volume,
                                epsilon );

        local_map->centroid( *cell_it, centroid() );
        TEST_FLOATING_EQUALITY( centroid[0], cell_width * i + cell_width / 2,
                                epsilon );
        TEST_FLOATING_EQUALITY( centroid[1], cell_width * j + cell_width / 2,
                                epsilon );
        TEST_FLOATING_EQUALITY( centroid[2], cell_width * k + cell_width / 2,
                                epsilon );

        TEST_ASSERT(
            local_map->isSafeToMapToReferenceFrame( *cell_it, centroid() ) );
        TEST_ASSERT(
            !local_map->isSafeToMapToReferenceFrame( *cell_it, bad_point() ) );

        TEST_ASSERT( local_map->mapToReferenceFrame( *cell_it, centroid(),
                                                     ref_point() ) );
        TEST_EQUALITY( ref_point[0], 0.0 );
        TEST_EQUALITY( ref_point[1], 0.0 );
        TEST_EQUALITY( ref_point[2], 0.0 );

        TEST_ASSERT( local_map->checkPointInclusion( *cell_it, ref_point() ) );
        TEST_ASSERT(
            !local_map->checkPointInclusion( *cell_it, bad_ref_point() ) );

        local_map->mapToPhysicalFrame( *cell_it, ref_point(), phys_point() );
        TEST_EQUALITY( phys_point[0], centroid[0] );
        TEST_EQUALITY( phys_point[1], centroid[1] );
        TEST_EQUALITY( phys_point[2], centroid[2] );

        shape_function->entitySupportIds( *cell_it, support_ids );
        TEST_EQUALITY( support_ids.size(), 8 );
        TEST_EQUALITY( support_ids[0], ( i ) + (j)*x_num_nodes +
                                           (k)*x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[1], ( i + 1 ) + (j)*x_num_nodes +
                                           (k)*x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[2], ( i + 1 ) + ( j + 1 ) * x_num_nodes +
                                           (k)*x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[3], ( i ) + ( j + 1 ) * x_num_nodes +
                                           (k)*x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[4],
                       ( i ) + (j)*x_num_nodes +
                           ( k + 1 ) * x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[5],
                       ( i + 1 ) + (j)*x_num_nodes +
                           ( k + 1 ) * x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[6],
                       ( i + 1 ) + ( j + 1 ) * x_num_nodes +
                           ( k + 1 ) * x_num_nodes * y_num_nodes );
        TEST_EQUALITY( support_ids[7],
                       ( i ) + ( j + 1 ) * x_num_nodes +
                           ( k + 1 ) * x_num_nodes * y_num_nodes );
    }

    // Check the nodal field.
    auto field = mesh.nodalField( 3 );
    auto field_supports = field->getLocalSupportIds();
    int counter = 0;
    for ( node_it = nodes_begin; node_it != nodes_end; ++node_it, ++counter )
    {
        TEST_EQUALITY( node_it->id(), field_supports[counter] );
    }
    for ( node_it = nodes_begin; node_it != nodes_end; ++node_it )
    {
        field->writeFieldData( node_it->id(), 0, node_it->id() + 32.2 );
        field->writeFieldData( node_it->id(), 1, node_it->id() + 3.2 );
    }

    for ( node_it = nodes_begin; node_it != nodes_end; ++node_it )
    {
        TEST_EQUALITY( field->readFieldData( node_it->id(), 0 ),
                       node_it->id() + 32.2 );
        TEST_EQUALITY( field->readFieldData( node_it->id(), 1 ),
                       node_it->id() + 3.2 );
    }
}

//---------------------------------------------------------------------------//
#define TRANSFORM( a00, a01, a02, a11, a12, a22, a0, a1, a2, a )               \
    ( a00 ) * c[0] * c[0] + (a01)*c[0] * c[1] + (a02)*c[0] * c[2] +            \
        (a11)*c[1] * c[1] + (a12)*c[1] * c[2] + (a22)*c[2] * c[2] +            \
        (a0)*c[0] + (a1)*c[1] + (a2)*c[2] + ( a )
#define TRANSFORM0( a00, a01, a02, a11, a12, a22, a0, a1, a2, a )              \
    2 * (a00)*c[0] + (a01)*c[1] + (a02)*c[2] + ( a0 )
#define TRANSFORM1( a00, a01, a02, a11, a12, a22, a0, a1, a2, a )              \
    2 * (a11)*c[1] + (a01)*c[0] + (a12)*c[2] + ( a1 )
#define TRANSFORM2( a00, a01, a02, a11, a12, a22, a0, a1, a2, a )              \
    2 * (a22)*c[2] + (a02)*c[0] + (a12)*c[1] + ( a2 )
void transform( const Teuchos::Array<double> &c,
                Teuchos::Array<double> &coords_out )
{
    DTK_REQUIRE( c.size() == 3 );
    DTK_REQUIRE( coords_out.size() == 3 );

    coords_out[0] =
        TRANSFORM( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.4, 9.1, 5.2, 1.2 );
    coords_out[1] =
        TRANSFORM( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.3, 2.6, 2.4, 6.2 );
    coords_out[2] =
        TRANSFORM( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1, 6.0, 2.3, 1.9 );
}

double transform_jacobian( int i, int j, const Teuchos::Array<double> &c )
{
    DTK_REQUIRE( i >= 0 && i < 3 );
    DTK_REQUIRE( j >= 0 && j < 3 );
    DTK_REQUIRE( c.size() == 3 );

    switch ( i )
    {
    case 0:
        switch ( j )
        {
        case 0:
            return TRANSFORM0( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.4, 9.1, 5.2,
                               1.2 );
        case 1:
            return TRANSFORM1( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.4, 9.1, 5.2,
                               1.2 );
        case 2:
            return TRANSFORM2( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.4, 9.1, 5.2,
                               1.2 );
        }
    case 1:
        switch ( j )
        {
        case 0:
            return TRANSFORM0( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.3, 2.6, 2.4,
                               6.2 );
        case 1:
            return TRANSFORM1( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.3, 2.6, 2.4,
                               6.2 );
        case 2:
            return TRANSFORM2( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.3, 2.6, 2.4,
                               6.2 );
        }
    case 2:
        switch ( j )
        {
        case 0:
            return TRANSFORM0( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1, 6.0, 2.3,
                               1.9 );
        case 1:
            return TRANSFORM1( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1, 6.0, 2.3,
                               1.9 );
        case 2:
            return TRANSFORM2( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1, 6.0, 2.3,
                               1.9 );
        }
    }
    return 0.0;
}

double transform_determinant( const Teuchos::Array<double> &c )
{
    Teuchos::Array<Teuchos::Array<double>> J( 3 );
    for ( int i = 0; i < 3; i++ )
    {
        J[i].resize( 3 );
        for ( int j = 0; j < 3; j++ )
            J[i][j] = transform_jacobian( i, j, c );
    }

    double det = J[0][0] * J[1][1] * J[2][2] + J[0][1] * J[1][2] * J[2][0] +
                 J[0][2] * J[1][0] * J[2][1] - J[0][2] * J[1][1] * J[2][0] -
                 J[0][1] * J[1][0] * J[2][2] - J[0][0] * J[1][2] * J[2][1];

    return det;
}
#undef TRANSORM
#undef TRANSFORM0
#undef TRANSFORM1
#undef TRANSFORM2

TEUCHOS_UNIT_TEST( ReferenceHexMesh, jacobian_test )
{
    // Get the comm.
    Teuchos::RCP<const Teuchos::Comm<int>> comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    if ( comm_size > 1 )
        return;

    // Create an entity set.
    Teuchos::RCP<DataTransferKit::BasicEntitySet> entity_set =
        Teuchos::rcp( new DataTransferKit::BasicEntitySet( comm, 3 ) );

    // Create the nodes.
    Teuchos::Array<DataTransferKit::Entity> hex_nodes( 8 );
    int node_owner = comm_rank;
    for ( int node_id = 0; node_id < 8; node_id++ )
    {
        Teuchos::Array<double> old_coords( 3 ), new_coords( 3 );
        old_coords[0] =
            ( node_id == 1 || node_id == 2 || node_id == 5 || node_id == 6 )
                ? 1.0
                : -1.0;
        old_coords[1] =
            ( node_id == 2 || node_id == 3 || node_id == 6 || node_id == 7 )
                ? 1.0
                : -1.0;
        old_coords[2] =
            ( node_id == 4 || node_id == 5 || node_id == 6 || node_id == 7 )
                ? 1.0
                : -1.0;

        transform( old_coords, new_coords );

        hex_nodes[node_id] = DataTransferKit::UnitTest::ReferenceNode(
            node_id, node_owner, new_coords[0], new_coords[1], new_coords[2] );

        // Add it to the entity set.
        entity_set->addEntity( hex_nodes[node_id] );
    }

    // Create the elements.
    DataTransferKit::Entity hex =
        DataTransferKit::UnitTest::ReferenceHex( 0, comm_rank, hex_nodes );

    // Add the element to the entity set.
    entity_set->addEntity( hex );

    // Create the function space.
    Teuchos::RCP<DataTransferKit::EntityLocalMap> lm =
        Teuchos::rcp( new DataTransferKit::UnitTest::ReferenceHexLocalMap() );
    Teuchos::RCP<DataTransferKit::EntityShapeFunction> sf = Teuchos::rcp(
        new DataTransferKit::UnitTest::ReferenceHexShapeFunction() );
    Teuchos::RCP<DataTransferKit::EntityIntegrationRule> ir = Teuchos::rcp(
        new DataTransferKit::UnitTest::ReferenceHexIntegrationRule() );
    Teuchos::RCP<DataTransferKit::FunctionSpace> function_space = Teuchos::rcp(
        new DataTransferKit::FunctionSpace( entity_set, lm, sf, ir ) );

    DataTransferKit::Jacobian J( function_space );

    Teuchos::Array<double> ref_point( 3 );
    ref_point[0] = 0.2;
    ref_point[1] = 0.7;
    ref_point[2] = 0.3;

    Teuchos::Array<Teuchos::Array<double>> jacobian;
    J.jacobian( hex, ref_point, jacobian );
    for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
            TEST_FLOATING_EQUALITY( jacobian[i][j],
                                    transform_jacobian( i, j, ref_point ),
                                    epsilon );
    TEST_FLOATING_EQUALITY( J.jacobian_determinant( hex, ref_point ),
                            transform_determinant( ref_point ), epsilon );
}

//---------------------------------------------------------------------------//
// end tstReferenceNode.cpp
//---------------------------------------------------------------------------//
