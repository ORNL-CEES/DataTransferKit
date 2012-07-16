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
 * \file DTK_TopologyTools.cpp
 * \author Stuart R. Slattery
 * \brief TopologyTools definition.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "DTK_CellTopologyFactory.hpp"
#include "DTK_TopologyTools.hpp"
#include <DTK_Exception.hpp>

#include <Teuchos_ENull.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_CellTools.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*! 
 * \brief Get the number of linear nodes for a particular iMesh topology.
 */
int TopologyTools::numLinearNodes( const moab::EntityType element_topology )
{
    int num_nodes = 0;
    
    switch( element_topology )
    {
	case moab::MBVERTEX:
	    num_nodes = 1;
	    break;

	case moab::MBEDGE:
	    num_nodes = 2;
	    break;

	case moab::MBTRI:
	    num_nodes = 3;
	    break;

	case moab::MBQUAD:
	    num_nodes = 4;
	    break;

	case moab::MBTET:
	    num_nodes = 4;
	    break;

	case moab::MBHEX:
	    num_nodes = 8;
	    break;

	case moab::MBPYRAMID:
	    num_nodes = 5;
	    break;

	default:
	    testPrecondition( moab::MBEDGE    == element_topology ||
			      moab::MBTRI     == element_topology ||
			      moab::MBQUAD    == element_topology ||
			      moab::MBTET     == element_topology ||
			      moab::MBHEX     == element_topology ||
			      moab::MBPYRAMID == element_topology,
			      "Invalid mesh topology" );
    }

    return num_nodes;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Reorder a list of element nodes from MBCN ordering to Shards
 * ordering. 
 */
void TopologyTools::MBCN2Shards( std::vector<moab::EntityHandle>& element_nodes, 
				 const int num_nodes,
				 const moab::EntityType topology )
{
    std::vector<moab::EntityHandle> temp_nodes( num_nodes );

    switch( topology )
    {
	case moab::MBHEX:

	    switch( num_nodes )
	    {
		case 27:

		    for ( int n = 0; n < 20; ++n )
		    {
			temp_nodes[n] = element_nodes[n];
		    }
		    temp_nodes[20] = element_nodes[26];
		    temp_nodes[21] = element_nodes[24];
		    temp_nodes[22] = element_nodes[25];
		    temp_nodes[23] = element_nodes[23];
		    temp_nodes[24] = element_nodes[21];
		    temp_nodes[25] = element_nodes[20];
		    temp_nodes[26] = element_nodes[22];

		    for ( int n = 0; n < num_nodes; ++n )
		    {
			element_nodes[n] = temp_nodes[n];
		    }
		    break;

		default: ;
	    }
	    break;

	default: ;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Point in element query.
 */
bool TopologyTools::pointInElement( Teuchos::Array<double>& coords,
				    const moab::EntityHandle element,
				    const Teuchos::RCP<moab::Interface>& moab )
{
    moab::ErrorCode error;

    // Wrap the point in a field container.
    int node_dim = coords.size();
    Teuchos::Tuple<int,2> point_dimensions;
    point_dimensions[0] = 1;
    point_dimensions[1] = node_dim;
    Teuchos::ArrayRCP<double> coords_view = Teuchos::arcpFromArray( coords );
    Intrepid::FieldContainer<double> point(
	Teuchos::Array<int>(point_dimensions), coords_view );

    // Get the element topology.
    moab::EntityType element_topology = moab->type_from_handle( element );

    // Get the element nodes.
    std::vector<moab::EntityHandle> element_nodes;
    error = moab->get_adjacencies( &element,
				   1,
				   0,
				   false,
				   element_nodes );
    testInvariant( moab::MB_SUCCESS == error, "Failure getting element nodes" );

    // Permute MBCN ordering to Shards ordering.
    int num_element_nodes = element_nodes.size();
    MBCN2Shards( element_nodes, num_element_nodes, element_topology );

    // Extract the node coordinates.
    Teuchos::Array<double> cell_node_coords( 3 * num_element_nodes );
    error = moab->get_coords( &element_nodes[0], 
			      element_nodes.size(), 
			      &cell_node_coords[0] );
    testInvariant( moab::MB_SUCCESS == error, 
		   "Failure getting node coordinates" );

    // Typical topology case.
    if ( moab::MBPYRAMID != element_topology )
    {
	// Create the Shards topology for the element type.
	Teuchos::RCP<shards::CellTopology> cell_topo = 
	    CellTopologyFactory::create( element_topology, num_element_nodes );

	// Reduce the dimension of the coordinates if necessary and wrap in a
	// field container. This means (for now at least) that 2D meshes must be
	// constructed from 2D nodes (this obviously won't work for 2D meshes that
	// have curvature).
	Teuchos::Tuple<int,3> cell_node_dimensions;
	cell_node_dimensions[0] = 1;
	cell_node_dimensions[1] = num_element_nodes;
	cell_node_dimensions[2] = node_dim;
	for ( int i = 2; i != node_dim-1 ; --i )
	{
	    for ( int n = element_nodes.size() - 1; n > -1; --n )
	    {
		cell_node_coords.erase( cell_node_coords.begin() + 3*n + i );
	    }
	}
	Intrepid::FieldContainer<double> cell_nodes( 
	    Teuchos::Array<int>(cell_node_dimensions), 
	    Teuchos::arcpFromArray( cell_node_coords ) );

	// Map the point to the reference frame of the cell.
	Intrepid::FieldContainer<double> reference_point( 1, node_dim );
	Intrepid::CellTools<double>::mapToReferenceFrame( reference_point,
							  point,
							  cell_nodes,
							  *cell_topo,
							  0 );

	// Check for reference point inclusion in the reference cell.
	return Intrepid::CellTools<double>::checkPointsetInclusion( 
	    reference_point, *cell_topo);
    }

    // We have to handle pyramids differently because Intrepid doesn't
    // currently support them with basis functions. Instead we'll resolve them
    // with two linear tetrahedrons and check for point inclusion in that set
    // instead.
    else
    {
	// Create the Shards topology for the linear tetrahedrons.
	Teuchos::RCP<shards::CellTopology> cell_topo = 
	    CellTopologyFactory::create( moab::MBTET, 4 );

	// Build 2 tetrahedrons from the 1 pyramid.
	testInvariant( node_dim == 3, "Pyramid elements must be 3D." );
	Intrepid::FieldContainer<double> cell_nodes( 1, 4, node_dim );

	// Tetrahederon 1.
	cell_nodes( 0, 0, 0 ) = cell_node_coords[0];
	cell_nodes( 0, 0, 1 ) = cell_node_coords[1];
	cell_nodes( 0, 0, 2 ) = cell_node_coords[2];

	cell_nodes( 0, 1, 0 ) = cell_node_coords[3];
	cell_nodes( 0, 1, 1 ) = cell_node_coords[4];
	cell_nodes( 0, 1, 2 ) = cell_node_coords[5];

	cell_nodes( 0, 2, 0 ) = cell_node_coords[6];
	cell_nodes( 0, 2, 1 ) = cell_node_coords[7];
	cell_nodes( 0, 2, 2 ) = cell_node_coords[8];

	cell_nodes( 0, 3, 0 ) = cell_node_coords[12];
	cell_nodes( 0, 3, 1 ) = cell_node_coords[13];
	cell_nodes( 0, 3, 2 ) = cell_node_coords[14];

	// Map the point to the reference frame of the first linear
	// tetrahedron. 
	Intrepid::FieldContainer<double> reference_point( 1, node_dim );
	Intrepid::CellTools<double>::mapToReferenceFrame( reference_point,
							  point,
							  cell_nodes,
							  *cell_topo,
							  0 );

	// Check for reference point inclusion in tetrahedron 1.
	bool in_tet_1 = Intrepid::CellTools<double>::checkPointsetInclusion( 
	    reference_point, *cell_topo);

	// If the point is the first tetrahedron, it is in the pyramid.
	if ( in_tet_1 )
	{
	    return true;
	}

	// Tetrahederon 2.
	cell_nodes( 0, 0, 0 ) = cell_node_coords[0];
	cell_nodes( 0, 0, 1 ) = cell_node_coords[1];
	cell_nodes( 0, 0, 2 ) = cell_node_coords[2];

	cell_nodes( 0, 1, 0 ) = cell_node_coords[6];
	cell_nodes( 0, 1, 1 ) = cell_node_coords[7];
	cell_nodes( 0, 1, 2 ) = cell_node_coords[8];

	cell_nodes( 0, 2, 0 ) = cell_node_coords[9];
	cell_nodes( 0, 2, 1 ) = cell_node_coords[10];
	cell_nodes( 0, 2, 2 ) = cell_node_coords[11];

	cell_nodes( 0, 3, 0 ) = cell_node_coords[12];
	cell_nodes( 0, 3, 1 ) = cell_node_coords[13];
	cell_nodes( 0, 3, 2 ) = cell_node_coords[14];

	// Map the point to the reference frame of the first linear
	// tetrahedron. 
	Intrepid::CellTools<double>::mapToReferenceFrame( reference_point,
							  point,
							  cell_nodes,
							  *cell_topo,
							  0 );

	// Check for reference point inclusion in tetrahedron 2.
	bool in_tet_2 = Intrepid::CellTools<double>::checkPointsetInclusion( 
	    reference_point, *cell_topo);

	// If the point is the second tetrahedron, it is in the pyramid.
	if ( in_tet_2 )
	{
	    return true;
	}
	else
	{
	    return false;
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_TopologyTools.cpp
//---------------------------------------------------------------------------//
