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

#include "reference_implementation/DTK_ReferenceHexMesh.hpp"

#include <DTK_BasicEntityPredicates.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>

#include <limits>

//---------------------------------------------------------------------------//
// TEST EPSILON
//---------------------------------------------------------------------------//

const double epsilon = 100.0 * std::numeric_limits<double>::epsilon();

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( ReferenceHexMesh, mesh_test )
{
    // Get the comm.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();
    
    // Create the mesh.
    int x_num_cells = 20;
    int y_num_cells = 12;
    int local_z_num_cells = 15;
    int z_num_cells = local_z_num_cells * comm_size;

    int local_num_cells = x_num_cells*y_num_cells*local_z_num_cells;
    int global_num_cells = x_num_cells*y_num_cells*z_num_cells;    

    int x_num_nodes = x_num_cells + 1;
    int y_num_nodes = y_num_cells + 1;
    int local_z_num_nodes = local_z_num_cells + 1;
    int z_num_nodes = z_num_cells + 1;

    int local_num_nodes = x_num_nodes*y_num_nodes;
    if ( comm_rank < comm_size - 1 )
    {
        local_num_nodes *= local_z_num_nodes - 1;
    }
    else
    {
        local_num_nodes *= local_z_num_nodes;
    }
    
    int global_num_nodes = x_num_nodes*y_num_nodes*z_num_nodes;

    Teuchos::Array<double> x_edges( x_num_nodes );
    Teuchos::Array<double> y_edges( y_num_nodes );
    Teuchos::Array<double> z_edges( z_num_nodes );

    double cell_width = 0.1;
    for ( int i = 0; i < x_num_nodes; ++i ) x_edges[i] = i*cell_width;
    for ( int i = 0; i < y_num_nodes; ++i ) y_edges[i] = i*cell_width;
    for ( int i = 0; i < z_num_nodes; ++i ) z_edges[i] = i*cell_width;    

    DataTransferKit::UnitTest::ReferenceHexMesh mesh(
        comm, x_edges, y_edges, z_edges );

    // Test the entity set.
    auto entity_set = mesh.functionSpace()->entitySet();
    TEST_EQUALITY( entity_set->physicalDimension(), 3 );

    Teuchos::Tuple<double,6> box;
    entity_set->localBoundingBox( box );
    TEST_EQUALITY( box[0], 0.0 );
    TEST_EQUALITY( box[1], 0.0 );
    TEST_EQUALITY( box[2], cell_width*local_z_num_cells*comm_rank );
    TEST_EQUALITY( box[3], cell_width*x_num_cells );
    TEST_EQUALITY( box[4], cell_width*y_num_cells );
    TEST_EQUALITY( box[5], cell_width*local_z_num_cells*(comm_rank+1) );

    entity_set->globalBoundingBox( box );
    TEST_EQUALITY( box[0], 0.0 );
    TEST_EQUALITY( box[1], 0.0 );
    TEST_EQUALITY( box[2], 0.0 );
    TEST_EQUALITY( box[3], cell_width*x_num_cells );
    TEST_EQUALITY( box[4], cell_width*y_num_cells );
    TEST_EQUALITY( box[5], cell_width*z_num_cells );

    // Check the nodes.
    DataTransferKit::Entity entity;
    for ( int k = local_z_num_cells*comm_rank;
          k < local_z_num_cells*comm_rank + local_z_num_nodes;
          ++k )
    {
        for ( int j = 0; j < y_num_nodes; ++j )
        {
            for ( int i = 0; i < x_num_nodes; ++i )
            {
                DataTransferKit::EntityId node_id =
                    i + j*x_num_nodes + k*x_num_nodes*y_num_nodes;
                TEST_ASSERT( node_id < global_num_nodes );
                
                entity_set->getEntity( node_id, 0, entity );

                entity.boundingBox( box );
                TEST_EQUALITY( box[0], i*cell_width );
                TEST_EQUALITY( box[1], j*cell_width );
                TEST_EQUALITY( box[2], k*cell_width );

                if ( k < local_z_num_cells*comm_rank + local_z_num_nodes - 1 ||
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
    for ( int k = local_z_num_cells*comm_rank;
          k < local_z_num_cells*(comm_rank+1);
          ++k )
    {
        for ( int j = 0; j < y_num_cells; ++j )
        {
            for ( int i = 0; i < x_num_cells; ++i )
            {
                DataTransferKit::EntityId cell_id =
                    i + j*x_num_cells + k*x_num_cells*y_num_cells;
                TEST_ASSERT( cell_id < global_num_cells );

                entity_set->getEntity( cell_id, 3, entity );

                entity.boundingBox( box );
                TEST_EQUALITY( box[0], i*cell_width );
                TEST_EQUALITY( box[1], j*cell_width );
                TEST_EQUALITY( box[2], k*cell_width );
                TEST_EQUALITY( box[3], (i+1)*cell_width );
                TEST_EQUALITY( box[4], (j+1)*cell_width );
                TEST_EQUALITY( box[5], (k+1)*cell_width );

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
        int k = node_it->id() / (x_num_nodes*y_num_nodes);
        int j = (node_it->id() - k*x_num_nodes*y_num_nodes) / x_num_nodes;
        int i = node_it->id() - j*x_num_nodes - k*x_num_nodes*y_num_nodes;
        TEST_EQUALITY( node_it->id(),
                       i + j*x_num_nodes + k*x_num_nodes*y_num_nodes );
        
        TEST_EQUALITY( local_map->measure(*node_it), 0.0 );

        local_map->centroid( *node_it, node_coords() );
        TEST_EQUALITY( node_coords[0], cell_width*i );
        TEST_EQUALITY( node_coords[1], cell_width*j );
        TEST_EQUALITY( node_coords[2], cell_width*k );

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
    double volume = cell_width*cell_width*cell_width;
    auto cells_begin = cell_it.begin();
    auto cells_end = cell_it.end();
    for ( cell_it = cells_begin; cell_it != cells_end; ++cell_it )
    {
        int k = cell_it->id() / (x_num_cells*y_num_cells);
        int j = (cell_it->id() - k*x_num_cells*y_num_cells) / x_num_cells;
        int i = cell_it->id() - j*x_num_cells - k*x_num_cells*y_num_cells;
        TEST_EQUALITY( cell_it->id(),
                       i + j*x_num_cells + k*x_num_cells*y_num_cells );

        TEST_FLOATING_EQUALITY( local_map->measure(*cell_it), volume, epsilon );
                
        local_map->centroid( *cell_it, centroid() );
        TEST_FLOATING_EQUALITY( centroid[0], cell_width*i+cell_width/2, epsilon );
        TEST_FLOATING_EQUALITY( centroid[1], cell_width*j+cell_width/2, epsilon );
        TEST_FLOATING_EQUALITY( centroid[2], cell_width*k+cell_width/2, epsilon );

        TEST_ASSERT( local_map->isSafeToMapToReferenceFrame(
                         *cell_it, centroid() ) );
        TEST_ASSERT( !local_map->isSafeToMapToReferenceFrame(
                         *cell_it, bad_point() ) );

        TEST_ASSERT( local_map->mapToReferenceFrame(
                         *cell_it, centroid(), ref_point() ) );
        TEST_EQUALITY( ref_point[0], 0.0 );
        TEST_EQUALITY( ref_point[1], 0.0 );
        TEST_EQUALITY( ref_point[2], 0.0 );

        TEST_ASSERT( local_map->checkPointInclusion(
                         *cell_it, ref_point() ) );
        TEST_ASSERT( !local_map->checkPointInclusion(
                         *cell_it, bad_ref_point() ) );

        local_map->mapToPhysicalFrame(
            *cell_it, ref_point(), phys_point() );
        TEST_EQUALITY( phys_point[0], centroid[0] );
        TEST_EQUALITY( phys_point[1], centroid[1] );
        TEST_EQUALITY( phys_point[2], centroid[2] );

        shape_function->entitySupportIds( *cell_it, support_ids );
        TEST_EQUALITY( support_ids.size(), 8 );
        TEST_EQUALITY(
            support_ids[0],
            (i) + (j)*x_num_nodes + (k)*x_num_nodes*y_num_nodes );
        TEST_EQUALITY(
            support_ids[1],
            (i+1) + (j)*x_num_nodes + (k)*x_num_nodes*y_num_nodes );
        TEST_EQUALITY(
            support_ids[2],
            (i+1) + (j+1)*x_num_nodes + (k)*x_num_nodes*y_num_nodes );
        TEST_EQUALITY(
            support_ids[3],
            (i) + (j+1)*x_num_nodes + (k)*x_num_nodes*y_num_nodes );
        TEST_EQUALITY(
            support_ids[4],
            (i) + (j)*x_num_nodes + (k+1)*x_num_nodes*y_num_nodes );
        TEST_EQUALITY(
            support_ids[5],
            (i+1) + (j)*x_num_nodes + (k+1)*x_num_nodes*y_num_nodes );
        TEST_EQUALITY(
            support_ids[6],
            (i+1) + (j+1)*x_num_nodes + (k+1)*x_num_nodes*y_num_nodes );
        TEST_EQUALITY(
            support_ids[7],
            (i) + (j+1)*x_num_nodes + (k+1)*x_num_nodes*y_num_nodes );
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
        field->writeFieldData( node_it->id(),
                               0,
                               node_it->id() + 32.2 );
        field->writeFieldData( node_it->id(),
                               1,
                               node_it->id() + 3.2 );
    }

    for ( node_it = nodes_begin; node_it != nodes_end; ++node_it )
    {
        TEST_EQUALITY( field->readFieldData(node_it->id(),0),
                       node_it->id() + 32.2 );
        TEST_EQUALITY( field->readFieldData(node_it->id(),1),
                       node_it->id() + 3.2 );
    }    
}

//---------------------------------------------------------------------------//
// end tstReferenceNode.cpp
//---------------------------------------------------------------------------//

