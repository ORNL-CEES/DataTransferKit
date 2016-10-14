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
 * \brief DTK_ReferenceHexMesh.cpp
 * \author Stuart R. Slattery
 * \brief Reference implementation of a parallel hex mesh for unit testing.
 */
//---------------------------------------------------------------------------//

#include "DTK_ReferenceHexMesh.hpp"
#include "DTK_ReferenceHex.hpp"
#include "DTK_ReferenceHexIntegrationRule.hpp"
#include "DTK_ReferenceHexLocalMap.hpp"
#include "DTK_ReferenceHexShapeFunction.hpp"
#include "DTK_ReferenceNode.hpp"

#include <DTK_BasicEntityPredicates.hpp>

#include <DTK_BasicEntitySet.hpp>
#include <DTK_EntityCenteredField.hpp>

#include <DTK_DBC.hpp>

namespace DataTransferKit
{
namespace UnitTest
{
//---------------------------------------------------------------------------//
// Num cells constructor.
ReferenceHexMesh::ReferenceHexMesh(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm, double x_min,
    double x_max, const int x_num_cells, double y_min, double y_max,
    const int y_num_cells, double z_min, double z_max, const int z_num_cells )
{
    Teuchos::Array<double> x_edges =
        buildEdgeArray( x_min, x_max, x_num_cells );
    Teuchos::Array<double> y_edges =
        buildEdgeArray( y_min, y_max, y_num_cells );
    Teuchos::Array<double> z_edges =
        buildEdgeArray( z_min, z_max, z_num_cells );
    buildMesh( comm, x_edges, y_edges, z_edges );
}

//---------------------------------------------------------------------------//
// Edge array constructor.
ReferenceHexMesh::ReferenceHexMesh(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
    const Teuchos::Array<double> &x_edges,
    const Teuchos::Array<double> &y_edges,
    const Teuchos::Array<double> &z_edges )
{
    buildMesh( comm, x_edges, y_edges, z_edges );
}

//---------------------------------------------------------------------------//
// Get the function space.
Teuchos::RCP<DataTransferKit::FunctionSpace>
ReferenceHexMesh::functionSpace() const
{
    return d_function_space;
}

//---------------------------------------------------------------------------//
// Create a nodal field over the locally-owned nodes.
Teuchos::RCP<DataTransferKit::Field>
ReferenceHexMesh::nodalField( const int field_dim ) const
{
    return createNodalField( field_dim, false );
}

//---------------------------------------------------------------------------//
// Create a nodal field over the locally-owned + ghosted nodes.
Teuchos::RCP<DataTransferKit::Field>
ReferenceHexMesh::ghostedNodalField( const int field_dim ) const
{
    return createNodalField( field_dim, true );
}

//---------------------------------------------------------------------------//
// Build the mesh.
void ReferenceHexMesh::buildMesh(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
    const Teuchos::Array<double> &x_edges,
    const Teuchos::Array<double> &y_edges,
    const Teuchos::Array<double> &z_edges )
{
    // Get comm parameters.
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Create the partitioning. Z partitioning only.
    DTK_REQUIRE( x_edges.size() > 1 );
    DTK_REQUIRE( y_edges.size() > 1 );
    DTK_REQUIRE( z_edges.size() > comm_size );

    int x_num_nodes = x_edges.size();
    int y_num_nodes = y_edges.size();
    int z_num_nodes = z_edges.size();
    DTK_REMEMBER( int total_nodes = x_num_nodes * y_num_nodes * z_num_nodes );

    int x_num_cells = x_num_nodes - 1;
    int y_num_cells = y_num_nodes - 1;
    int z_num_cells = z_num_nodes - 1;
    DTK_REMEMBER( int total_cells = x_num_cells * y_num_cells * z_num_cells );

    Teuchos::Array<int> local_z_num_cells( comm_size, z_num_cells / comm_size );
    DTK_REMEMBER( int total_z = 0 );
    for ( int p = 0; p < comm_size; ++p )
    {
        if ( p < ( z_num_cells % comm_size ) )
            ++local_z_num_cells[p];
        DTK_REMEMBER( total_z += local_z_num_cells[p] );
    }
    DTK_CHECK( z_num_cells == total_z );

    Teuchos::Array<int> z_offsets( comm_size, 0 );
    for ( int p = 1; p < comm_size; ++p )
    {
        z_offsets[p] += local_z_num_cells[p - 1] + z_offsets[p - 1];
    }

    // Create an entity set.
    Teuchos::RCP<DataTransferKit::BasicEntitySet> entity_set =
        Teuchos::rcp( new DataTransferKit::BasicEntitySet( comm, 3 ) );

    // Create the nodes.
    int node_id = 0;
    int node_owner = 0;
    for ( int k = z_offsets[comm_rank];
          k < z_offsets[comm_rank] + local_z_num_cells[comm_rank] + 1; ++k )
    {
        for ( int j = 0; j < y_num_nodes; ++j )
        {
            for ( int i = 0; i < x_num_nodes; ++i )
            {
                // Create the node id.
                node_id = i + j * x_num_nodes + k * x_num_nodes * y_num_nodes;
                DTK_CHECK( node_id < total_nodes );

                // Get the owner rank of the node. This comm rank has
                // ownership of all nodes that compose local cells that do not
                // lie on the high z boundary. The last comm rank does own the
                // high z boundary because it has no neighbor with greater z
                // values.
                if ( k < z_offsets[comm_rank] + local_z_num_cells[comm_rank] ||
                     comm_rank == comm_size - 1 )
                {
                    node_owner = comm_rank;
                }
                else
                {
                    node_owner = comm_rank + 1;
                }

                // Create the node.
                DataTransferKit::Entity node = ReferenceNode(
                    node_id, node_owner, x_edges[i], y_edges[j], z_edges[k] );

                // Add it to the entity set.
                entity_set->addEntity( node );
            }
        }
    }

    // Create the elements.
    Teuchos::Array<DataTransferKit::Entity> hex_nodes( 8 );
    int element_id = 0;
    for ( int k = z_offsets[comm_rank];
          k < z_offsets[comm_rank] + local_z_num_cells[comm_rank]; ++k )
    {
        for ( int j = 0; j < y_num_cells; ++j )
        {
            for ( int i = 0; i < x_num_cells; ++i )
            {
                // Create the element id.
                element_id =
                    i + j * x_num_cells + k * x_num_cells * y_num_cells;
                DTK_CHECK( element_id < total_cells );

                // node 0
                node_id =
                    ( i ) + (j)*x_num_nodes + (k)*x_num_nodes * y_num_nodes;
                DTK_CHECK( node_id < total_nodes );
                entity_set->getEntity( node_id, 0, hex_nodes[0] );

                // node 1
                node_id =
                    ( i + 1 ) + (j)*x_num_nodes + (k)*x_num_nodes * y_num_nodes;
                DTK_CHECK( node_id < total_nodes );
                entity_set->getEntity( node_id, 0, hex_nodes[1] );

                // node 2
                node_id = ( i + 1 ) + ( j + 1 ) * x_num_nodes +
                          (k)*x_num_nodes * y_num_nodes;
                DTK_CHECK( node_id < total_nodes );
                entity_set->getEntity( node_id, 0, hex_nodes[2] );

                // node 3
                node_id = ( i ) + ( j + 1 ) * x_num_nodes +
                          (k)*x_num_nodes * y_num_nodes;
                DTK_CHECK( node_id < total_nodes );
                entity_set->getEntity( node_id, 0, hex_nodes[3] );

                // node 4
                node_id = ( i ) + (j)*x_num_nodes +
                          ( k + 1 ) * x_num_nodes * y_num_nodes;
                DTK_CHECK( node_id < total_nodes );
                entity_set->getEntity( node_id, 0, hex_nodes[4] );

                // node 5
                node_id = ( i + 1 ) + (j)*x_num_nodes +
                          ( k + 1 ) * x_num_nodes * y_num_nodes;
                DTK_CHECK( node_id < total_nodes );
                entity_set->getEntity( node_id, 0, hex_nodes[5] );

                // node 6
                node_id = ( i + 1 ) + ( j + 1 ) * x_num_nodes +
                          ( k + 1 ) * x_num_nodes * y_num_nodes;
                DTK_CHECK( node_id < total_nodes );
                entity_set->getEntity( node_id, 0, hex_nodes[6] );

                // node 7
                node_id = ( i ) + ( j + 1 ) * x_num_nodes +
                          ( k + 1 ) * x_num_nodes * y_num_nodes;
                DTK_CHECK( node_id < total_nodes );
                entity_set->getEntity( node_id, 0, hex_nodes[7] );

                // Create the element.
                DataTransferKit::Entity hex =
                    ReferenceHex( element_id, comm_rank, hex_nodes );

                // Add the element to the entity set.
                entity_set->addEntity( hex );
            }
        }
    }

    // Create the function space.
    Teuchos::RCP<DataTransferKit::EntityLocalMap> lm =
        Teuchos::rcp( new ReferenceHexLocalMap() );
    Teuchos::RCP<DataTransferKit::EntityShapeFunction> sf =
        Teuchos::rcp( new ReferenceHexShapeFunction() );
    Teuchos::RCP<DataTransferKit::EntityIntegrationRule> ir =
        Teuchos::rcp( new ReferenceHexIntegrationRule() );
    d_function_space = Teuchos::rcp(
        new DataTransferKit::FunctionSpace( entity_set, lm, sf, ir ) );
}

//---------------------------------------------------------------------------//
// Build an edge array.
Teuchos::Array<double>
ReferenceHexMesh::buildEdgeArray( const double min, const double max,
                                  const int num_cells ) const
{
    int num_nodes = num_cells + 1;
    double cell_width = ( max - min ) / num_cells;
    Teuchos::Array<double> edges( num_nodes );
    for ( int n = 0; n < num_nodes; ++n )
    {
        edges[n] = min + n * cell_width;
    }
    DTK_CHECK( std::abs( edges.back() - max ) < 1.0e-6 );
    return edges;
}

//---------------------------------------------------------------------------//
// Create a nodal field over the locally-owned nodes.
Teuchos::RCP<DataTransferKit::Field>
ReferenceHexMesh::createNodalField( const int field_dim,
                                    const bool is_ghosted ) const
{
    // Local + ghosted predicate.
    DataTransferKit::SelectAllPredicate all_pred;

    // Local only predicate.
    DataTransferKit::LocalEntityPredicate local_pred(
        d_function_space->entitySet()->communicator()->getRank() );

    // Select the nodes.
    DataTransferKit::PredicateFunction select_function =
        ( is_ghosted ) ? all_pred.getFunction() : local_pred.getFunction();
    auto node_it =
        d_function_space->entitySet()->entityIterator( 0, select_function );

    // Extract the node ids.
    auto nodes_begin = node_it.begin();
    auto nodes_end = node_it.end();
    Teuchos::Array<DataTransferKit::EntityId> node_ids;
    node_ids.reserve( node_it.size() );
    for ( node_it = nodes_begin; node_it != nodes_end; ++node_it )
    {
        node_ids.push_back( node_it->id() );
    }

    // Build the field.
    Teuchos::ArrayRCP<double> field_data( node_it.size() * field_dim );
    return Teuchos::rcp( new DataTransferKit::EntityCenteredField(
        node_ids(), field_dim, field_data,
        DataTransferKit::EntityCenteredField::BLOCKED ) );
}

//---------------------------------------------------------------------------//

} // end namespace UnitTest
} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_ReferenceHexMesh.hpp
//---------------------------------------------------------------------------//
