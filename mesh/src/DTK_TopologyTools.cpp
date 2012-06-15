//---------------------------------------------------------------------------//
/*!
 * \file DTK_TopologyTools.cpp
 * \author Stuart R. Slattery
 * \brief TopologyTools definition
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
 * \brief Point in element query.
 */
bool TopologyTools::pointInElement( Teuchos::Array<double>& coords,
				    const moab::EntityHandle element,
				    const Teuchos::RCP<moab::Interface>& moab )
{
    moab::ErrorCode error;

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

    // Extract only the nodes to build the linear element. This will need to
    // be updated for the higher order topologies.
    int num_linear_nodes = numLinearNodes( element_topology );
    Teuchos::Array<moab::EntityHandle> linear_nodes;
    for ( int n = 0; n < num_linear_nodes; ++n )
    {
	linear_nodes.push_back( element_nodes[n] );
    }

    // Create the Shards topology for the element type.
    Teuchos::RCP<shards::CellTopology> cell_topo = 
	CellTopologyFactory::create( element_topology, num_linear_nodes );

    // Extract the node coordinates.
    Teuchos::Array<double> cell_node_coords( 3 * num_linear_nodes );
    error = moab->get_coords( &linear_nodes[0], 
			      linear_nodes.size(), 
			      &cell_node_coords[0] );
    testInvariant( moab::MB_SUCCESS == error, 
		   "Failure getting node coordinates" );

    // Reduce the dimension of the coordinates if necessary and wrap in a
    // field container. This means (for now at least) that 2D meshes must be
    // constructed from 2D nodes (this obviously won't work for 2D meshes that
    // have curvature) 
    int node_dim = coords.size();
    Teuchos::Tuple<int,3> cell_node_dimensions;
    cell_node_dimensions[0] = 1;
    cell_node_dimensions[1] = num_linear_nodes;
    cell_node_dimensions[2] = node_dim;
    for ( int i = 2; i != node_dim-1 ; --i )
    {
	for ( int n = linear_nodes.size() - 1; n > -1; --n )
	{
	    cell_node_coords.erase( cell_node_coords.begin() + 3*n + i );
	}
    }
    Intrepid::FieldContainer<double> cell_nodes( 
	Teuchos::Array<int>(cell_node_dimensions), 
	Teuchos::arcpFromArray( cell_node_coords ) );

    // Wrap the point in a field container.
    Teuchos::Tuple<int,2> point_dimensions;
    point_dimensions[0] = 1;
    point_dimensions[1] = node_dim;
    Teuchos::ArrayRCP<double> coords_view = Teuchos::arcpFromArray( coords );
    Intrepid::FieldContainer<double> point(
	Teuchos::Array<int>(point_dimensions), coords_view );

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

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_TopologyTools.cpp
//---------------------------------------------------------------------------//
