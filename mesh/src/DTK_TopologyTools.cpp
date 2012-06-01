//---------------------------------------------------------------------------//
/*!
 * \file TopologyTools.cpp
 * \author Stuart Slattery
 * \brief TopologyTools definition
 */
//---------------------------------------------------------------------------//

#include <vector>
#include <cassert>

#include "DTK_CellTopologyFactory.hpp"
#include "DTK_TopologyTools.hpp"
#include <DTK_Exception.hpp>

#include <Teuchos_ENull.hpp>
#include <Teuchos_RCP.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_CellTools.hpp>

namespace FOOD
{

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

	default:
	    testPrecondition( moab::MBEDGE == element_topology ||
			      moab::MBTRI  == element_topology ||
			      moab::MBQUAD == element_topology ||
			      moab::MBTET  == element_topology ||
			      moab::MBHEX  == element_topology ,
			      "Invalid mesh topology" );
    }

    return num_nodes;
}

/*!
 * \brief Point in volume query on an element
 */
bool TopologyTools::pointInVolume( const double coords[3],
				   const moab::EntityHandle element,
				   const Teuchos::RCP<moab::Interface>& moab )
{
    moab::ErrorCode error;
    int topology = 0;
    iMesh_getEntTopo( mesh, entity, &topology, &error );
    assert( iBase_SUCCESS == error );

    iBase_EntityHandle *element_nodes = 0;
    int element_nodes_allocated = 0;
    int element_nodes_size = 0;
    iMesh_getEntAdj( mesh,
		     entity,
		     iBase_VERTEX,
		     &element_nodes,
		     &element_nodes_allocated,
		     &element_nodes_size,
		     &error );
    assert( iBase_SUCCESS == error );

    MBCN2Shards( element_nodes, element_nodes_size, topology );

    int num_linear_nodes = numLinearNodes( topology );

    CellTopologyFactory topo_factory;
    Teuchos::RCP<shards::CellTopology> cell_topo = 
	topo_factory.create( topology, num_linear_nodes );

    std::vector<iBase_EntityHandle> linear_nodes;
    for ( int n = 0; n < num_linear_nodes; ++n )
    {
	linear_nodes.push_back( element_nodes[n] );
    }

    int coords_allocated = 0;
    int coords_size = 0;
    double *coord_array = 0;
    iMesh_getVtxArrCoords( mesh,
			   &linear_nodes[0],
			   num_linear_nodes,
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    Teuchos::Tuple<int,3> cell_node_dimensions;
    cell_node_dimensions[0] = 1;
    cell_node_dimensions[1] = num_linear_nodes;
    cell_node_dimensions[2] = 3;
    Intrepid::FieldContainer<double> cell_nodes( 
	Teuchos::Array<int>(cell_node_dimensions), coord_array );

    Intrepid::FieldContainer<double> find_coords(1,3);
    find_coords(0,0) = coords[0];
    find_coords(0,1) = coords[1];
    find_coords(0,2) = coords[2];
    Intrepid::FieldContainer<double> reference_points(1,3);
    Intrepid::CellTools<double>::mapToReferenceFrame( reference_points,
						      find_coords,
						      cell_nodes,
						      *cell_topo,
						      0 );

    bool return_val = Intrepid::CellTools<double>::checkPointInclusion( 
	&reference_points[0],
	3,
	*cell_topo);

    free( element_nodes );
    free( coord_array );
    cell_topo = Teuchos::null;

    return return_val;
}

} // end namespace FOOD

//---------------------------------------------------------------------------//
// end TopologyTools.cpp
//---------------------------------------------------------------------------//
