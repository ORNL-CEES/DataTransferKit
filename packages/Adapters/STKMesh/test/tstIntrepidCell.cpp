//---------------------------------------------------------------------------//
/*!
 * \file   tstIntrepidCell
 * \author Stuart Slattery
 * \brief  IntrepidCell and IntrepidSideCell tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <sstream>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Comm.hpp>

#include <Shards_CellTopology.hpp>
#include <Shards_BasicTopologies.hpp>

#include <Intrepid_FieldContainer.hpp>

#include "DTK_IntrepidCell.hpp"
#include "DTK_IntrepidSideCell.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
// Get the default communicator.
template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
    return Teuchos::DefaultComm<Ordinal>::getComm();
}

//---------------------------------------------------------------------------//
// Given an element id compute the ids of the associated nodes.
void compute_elem_node_ids( int elem_id, int node_ids[] )
{
    const unsigned nodes_per_side = 4;
    const unsigned base = elem_id * nodes_per_side ;
    node_ids[0] = base ;
    node_ids[1] = base + 4 ;
    node_ids[2] = base + 5 ;
    node_ids[3] = base + 1 ;
    node_ids[4] = base + 3 ;
    node_ids[5] = base + 7 ;
    node_ids[6] = base + 6 ;
    node_ids[7] = base + 2 ;
}

//---------------------------------------------------------------------------//
// Given a node_id compute its spatial coordinates
void compute_node_coordinates( int node_id, double coord[] )
{
    // i_length is the same as the number of the side it occurs in
    // i_plane is the position of the node in the side
    const unsigned i_length = node_id / 4 ;
    const unsigned i_plane  = node_id % 4 ;

    coord[0] = i_length ;
    coord[1] = i_plane == 1 || i_plane == 2 ? 1.0 : 0.0 ;
    coord[2] = i_plane == 2 || i_plane == 3 ? 1.0 : 0.0 ;
}

//---------------------------------------------------------------------------//
// Build the mesh node coordinates.
Teuchos::Array<double> buildMeshNodeCoordinates(
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const int dimension, const int num_elements )
{
    // Generate the node ids.
    std::set<int> node_set;
    int dim = 3;
    const int nodes_per_elem = shards::Hexahedron<8>::node_count;
    int node_ids[nodes_per_elem];
    for ( int i = 0; i < num_elements; ++i )
    {
	compute_elem_node_ids( i, node_ids );
	for ( int j = 0; j < nodes_per_elem; ++j )
	{
	    node_set.insert( node_ids[j] );
	}
    }

    // Generate the node coordinates.
    int num_nodes = node_set.size();
    Teuchos::Array<double> coord_array( dim * num_nodes );   
    for ( int i = 0; i < num_nodes; ++i )
    {
	compute_node_coordinates( i, coord_array(dim*i,dim).getRawPtr() );
    }

    // Create the unrolled element node coordinates.
    Teuchos::Array<double> elem_coord_array( 
	dim*num_elements*nodes_per_elem );
    for ( int i = 0; i < num_elements; ++i )
    {
	compute_elem_node_ids( i, node_ids );
	for ( int j = 0; j < nodes_per_elem; ++j )
	{
	    for ( int d = 0; d < dim; ++d )
	    {
		elem_coord_array[ i*nodes_per_elem*dim + j*dim + d ]
		    = coord_array[ node_ids[j]*dim + d ];
	    }
	}
    }

    return elem_coord_array;
}

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( IntrepidCell, element_all_cells_test )
{
    // INITIALIZATION
    // --------------

    // Basic problem info.
    int dimension = 3;
    int num_elements = 5;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();

    // Build the mesh.
    Teuchos::Array<double> element_coordinates =
	buildMeshNodeCoordinates( comm, dimension, num_elements );

    // Build a cell manager for the mesh.
    shards::CellTopology element_topo =
	shards::getCellTopologyData<shards::Hexahedron<8> >();
    int degree = 1;
    DataTransferKit::IntrepidCell intrepid_cell( element_topo, degree );

    // For each element in the mesh compute cell measures, integrals, and
    // physical cubature points and check them.
    int workset_size = num_elements;
    int num_ip = 1;
    int num_nodes = element_topo.getNodeCount();
    Teuchos::Array<int> coord_dims( 3 );
    coord_dims[0] = workset_size;
    coord_dims[1] = num_nodes;
    coord_dims[2] = dimension;
    Intrepid::FieldContainer<double> element_coords( 
	coord_dims, element_coordinates() );
    Intrepid::FieldContainer<double> cell_measure(workset_size);
    Intrepid::FieldContainer<double> ip_coords;
    Intrepid::FieldContainer<double> integrals(workset_size);
    Intrepid::FieldContainer<double> dofs(workset_size,num_ip);
    int dof_val = 2.0;
    dofs.initialize( dof_val );
    Intrepid::FieldContainer<double> param_coords( 1, dimension );
    param_coords.initialize(0.0);
    Intrepid::FieldContainer<double> physical_coords( workset_size, 1, dimension );

    intrepid_cell.allocateCellState( element_coords );
    intrepid_cell.updateCellState();

    intrepid_cell.mapToCellPhysicalFrame( param_coords, physical_coords );
    intrepid_cell.getCellMeasures( cell_measure );
    intrepid_cell.getPhysicalIntegrationCoordinates( ip_coords );
    intrepid_cell.integrate( dofs, integrals );

    TEST_EQUALITY( workset_size, intrepid_cell.getNumCells() );
    TEST_EQUALITY( num_ip, intrepid_cell.getNumIntegrationPoints() );
    TEST_EQUALITY( dimension, intrepid_cell.getSpatialDimension() );

    for ( int cell = 0; cell < num_elements; ++cell )
    {
	for ( int ip = 0; ip < ip_coords.dimension(1); ++ip )
	{
	    TEST_EQUALITY( ip_coords(cell,ip,0), 1.0*(cell) + 0.5 );
	    TEST_EQUALITY( ip_coords(cell,ip,1), 0.5 );
	    TEST_EQUALITY( ip_coords(cell,ip,2), 0.5 );
	}
	
	TEST_EQUALITY( physical_coords(cell,0,0), 1.0*(cell) + 0.5 );
	TEST_EQUALITY( physical_coords(cell,0,1), 0.5 );
	TEST_EQUALITY( physical_coords(cell,0,2), 0.5 );

	TEST_EQUALITY( dof_val*cell_measure(cell), integrals(cell) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( IntrepidCell, element_single_cell_test )
{
    // INITIALIZATION
    // --------------

    // Basic problem info.
    int dimension = 3;
    int num_elements = 5;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();

    // Build the mesh
    Teuchos::Array<double> element_coordinates =
	buildMeshNodeCoordinates( comm, dimension, num_elements );

    // Build an cell manager for the mesh.
    shards::CellTopology element_topo =
	shards::getCellTopologyData<shards::Hexahedron<8> >();
    int degree = 1;
    DataTransferKit::IntrepidCell intrepid_cell( element_topo, degree );

    // For each element in the mesh compute cell measures, integrals, and
    // physical cubature points and check them.
    int workset_size = 1;
    int num_ip = 1;
    int num_nodes = element_topo.getNodeCount();
    Teuchos::Array<int> coord_dims( 3 );
    coord_dims[0] = workset_size;
    coord_dims[1] = num_nodes;
    coord_dims[2] = dimension;
    Intrepid::FieldContainer<double> cell_measure(workset_size);
    Intrepid::FieldContainer<double> ip_coords;
    Intrepid::FieldContainer<double> integrals(workset_size);
    Intrepid::FieldContainer<double> dofs(workset_size,num_ip);
    int dof_val = 2.0;
    dofs.initialize( dof_val );
    Intrepid::FieldContainer<double> param_coords( 1, dimension );
    param_coords.initialize(0.0);
    Intrepid::FieldContainer<double> physical_coords( workset_size, 1, dimension );

    for ( int cell = 0; cell < num_elements; ++cell )
    {
	Intrepid::FieldContainer<double> element_coords( 
	    coord_dims, element_coordinates(num_nodes*cell*dimension,dimension*num_nodes) );
	intrepid_cell.allocateCellState( element_coords );
	intrepid_cell.updateCellState();

	TEST_EQUALITY( workset_size, intrepid_cell.getNumCells() );
	TEST_EQUALITY( num_ip, intrepid_cell.getNumIntegrationPoints() );
	TEST_EQUALITY( dimension, intrepid_cell.getSpatialDimension() );

	intrepid_cell.mapToCellPhysicalFrame( param_coords, physical_coords );
	TEST_EQUALITY( physical_coords(0,0,0), 1.0*(cell) + 0.5 );
	TEST_EQUALITY( physical_coords(0,0,1), 0.5 );
	TEST_EQUALITY( physical_coords(0,0,2), 0.5 );

	intrepid_cell.getCellMeasures( cell_measure );
	intrepid_cell.getPhysicalIntegrationCoordinates( ip_coords );
	for ( int ip = 0; ip < ip_coords.dimension(1); ++ip )
	{
	    TEST_EQUALITY( ip_coords(0,ip,0), 1.0*(cell) + 0.5 );
	    TEST_EQUALITY( ip_coords(0,ip,1), 0.5 );
	    TEST_EQUALITY( ip_coords(0,ip,2), 0.5 );
	}

	intrepid_cell.integrate( dofs, integrals );
	TEST_EQUALITY( dof_val*cell_measure(0), integrals(0) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( IntrepidCell, side_all_cells_test )
{
    // INITIALIZATION
    // --------------

    // Basic problem info.
    int dimension = 3;
    int num_elements = 5;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();

    // Build the mesh
    Teuchos::Array<double> element_coordinates =
	buildMeshNodeCoordinates( comm, dimension, num_elements );

    // Build an cell manager for the mesh.
    shards::CellTopology element_topo =
	shards::getCellTopologyData<shards::Hexahedron<8> >();
    shards::CellTopology side_topo =
	shards::getCellTopologyData<shards::Quadrilateral<4> >();
    int degree = 1;
    int side_id = 0;
    DataTransferKit::IntrepidSideCell 
	intrepid_cell( side_topo, side_id, element_topo, degree);

    // For each element in the mesh compute cell measures, integrals, and
    // physical cubature points and check them.
    int workset_size = num_elements;
    int num_ip = 1;
    int num_nodes = element_topo.getNodeCount();
    Teuchos::Array<int> coord_dims( 3 );
    coord_dims[0] = workset_size;
    coord_dims[1] = num_nodes;
    coord_dims[2] = dimension;
    Intrepid::FieldContainer<double> element_coords( 
	coord_dims, element_coordinates() );
    Intrepid::FieldContainer<double> cell_measure(workset_size);
    Intrepid::FieldContainer<double> ip_coords;
    Intrepid::FieldContainer<double> integrals(workset_size);
    Intrepid::FieldContainer<double> dofs(workset_size,num_ip);
    int dof_val = 2.0;
    dofs.initialize( dof_val );
    Intrepid::FieldContainer<double> param_coords( 1, dimension-1 );
    param_coords.initialize(0.0);
    Intrepid::FieldContainer<double> physical_coords( workset_size, 1, dimension );

    intrepid_cell.allocateCellState( element_coords );
    intrepid_cell.updateCellState();

    intrepid_cell.mapToCellPhysicalFrame( param_coords, physical_coords );
    intrepid_cell.getCellMeasures( cell_measure );
    intrepid_cell.getPhysicalIntegrationCoordinates( ip_coords );
    intrepid_cell.integrate( dofs, integrals );

    TEST_EQUALITY( workset_size, intrepid_cell.getNumCells() );
    TEST_EQUALITY( num_ip, intrepid_cell.getNumIntegrationPoints() );
    TEST_EQUALITY( dimension, intrepid_cell.getSpatialDimension() );

    for ( int cell = 0; cell < num_elements; ++cell )
    {
	for ( int ip = 0; ip < ip_coords.dimension(1); ++ip )
	{
	    TEST_EQUALITY( ip_coords(cell,ip,0), 1.0*(cell) + 0.5 );
	    TEST_EQUALITY( ip_coords(cell,ip,1), 0.0 );
	    TEST_EQUALITY( ip_coords(cell,ip,2), 0.5 );
	}

	TEST_EQUALITY( physical_coords(cell,0,0), 1.0*(cell) + 0.5 );
	TEST_EQUALITY( physical_coords(cell,0,1), 0.0 );
	TEST_EQUALITY( physical_coords(cell,0,2), 0.5 );

	TEST_EQUALITY( dof_val*cell_measure(cell), integrals(cell) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( IntrepidCell, side_single_cell_test )
{
    // INITIALIZATION
    // --------------

    // Basic problem info.
    int dimension = 3;
    int num_elements = 5;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();

    // Build the mesh
    Teuchos::Array<double> element_coordinates =
	buildMeshNodeCoordinates( comm, dimension, num_elements );

    // Build an cell manager for the mesh.
    shards::CellTopology element_topo =
	shards::getCellTopologyData<shards::Hexahedron<8> >();
    shards::CellTopology side_topo =
	shards::getCellTopologyData<shards::Quadrilateral<4> >();
    int degree = 1;
    int side_id = 0;
    DataTransferKit::IntrepidSideCell
	intrepid_cell( side_topo, side_id, element_topo, degree);

    // For each element in the mesh compute cell measures, integrals, and
    // physical cubature points and check them.
    int workset_size = 1;
    int num_ip = 1;
    int num_nodes = element_topo.getNodeCount();
    Teuchos::Array<int> coord_dims( 3 );
    coord_dims[0] = workset_size;
    coord_dims[1] = num_nodes;
    coord_dims[2] = dimension;
    Intrepid::FieldContainer<double> cell_measure(workset_size);
    Intrepid::FieldContainer<double> ip_coords;
    Intrepid::FieldContainer<double> integrals(workset_size);
    Intrepid::FieldContainer<double> dofs(workset_size,num_ip);
    int dof_val = 2.0;
    dofs.initialize( dof_val );
    Intrepid::FieldContainer<double> param_coords( 1, dimension-1 );
    param_coords.initialize(0.0);
    Intrepid::FieldContainer<double> physical_coords( workset_size, 1, dimension );

    for ( int cell = 0; cell < num_elements; ++cell )
    {
	Intrepid::FieldContainer<double> element_coords( 
	    coord_dims, element_coordinates(num_nodes*cell*dimension,dimension*num_nodes) );
	intrepid_cell.allocateCellState( element_coords );
	intrepid_cell.updateCellState();

	TEST_EQUALITY( workset_size, intrepid_cell.getNumCells() );
	TEST_EQUALITY( num_ip, intrepid_cell.getNumIntegrationPoints() );
	TEST_EQUALITY( dimension, intrepid_cell.getSpatialDimension() );

	intrepid_cell.mapToCellPhysicalFrame( param_coords, physical_coords );
	TEST_EQUALITY( physical_coords(0,0,0), 1.0*(cell) + 0.5 );
	TEST_EQUALITY( physical_coords(0,0,1), 0.0 );
	TEST_EQUALITY( physical_coords(0,0,2), 0.5 );

	intrepid_cell.getCellMeasures( cell_measure );
	intrepid_cell.getPhysicalIntegrationCoordinates( ip_coords );
	for ( int ip = 0; ip < ip_coords.dimension(1); ++ip )
	{
	    TEST_EQUALITY( ip_coords(0,ip,0), 1.0*(cell) + 0.5 );
	    TEST_EQUALITY( ip_coords(0,ip,1), 0.0 );
	    TEST_EQUALITY( ip_coords(0,ip,2), 0.5 );
	}

	intrepid_cell.integrate( dofs, integrals );
	TEST_EQUALITY( dof_val*cell_measure(0), integrals(0) );
    }
}
 
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( IntrepidCell, element_ffupdate_all_cells_test )
{
    // INITIALIZATION
    // --------------

    // Basic problem info.
    int dimension = 3;
    int num_elements = 5;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();

    // Build the mesh
    Teuchos::Array<double> element_coordinates =
	buildMeshNodeCoordinates( comm, dimension, num_elements );

    // Build an cell manager for the mesh.
    shards::CellTopology element_topo =
	shards::getCellTopologyData<shards::Hexahedron<8> >();
    int degree = 1;
    DataTransferKit::IntrepidCell intrepid_cell( element_topo, degree );

    // For each element in the mesh compute cell measures, integrals, and
    // physical cubature points and check them.
    int workset_size = num_elements;
    int num_ip = 1;
    int num_nodes = element_topo.getNodeCount();
    Teuchos::Array<int> coord_dims( 3 );
    coord_dims[0] = workset_size;
    coord_dims[1] = num_nodes;
    coord_dims[2] = dimension;
    Intrepid::FieldContainer<double> element_coords( 
	coord_dims, element_coordinates() );
    Intrepid::FieldContainer<double> cell_measure(workset_size);
    Intrepid::FieldContainer<double> ip_coords;
    Intrepid::FieldContainer<double> integrals(workset_size);
    Intrepid::FieldContainer<double> dofs(workset_size,num_ip);
    int dof_val = 2.0;
    dofs.initialize( dof_val );

    DataTransferKit::IntrepidCell::updateState(
	intrepid_cell, element_coords );

    intrepid_cell.getCellMeasures( cell_measure );
    intrepid_cell.getPhysicalIntegrationCoordinates( ip_coords );
    intrepid_cell.integrate( dofs, integrals );

    TEST_EQUALITY( workset_size, intrepid_cell.getNumCells() );
    TEST_EQUALITY( num_ip, intrepid_cell.getNumIntegrationPoints() );
    TEST_EQUALITY( dimension, intrepid_cell.getSpatialDimension() );

    for ( int cell = 0; cell < num_elements; ++cell )
    {
	for ( int ip = 0; ip < ip_coords.dimension(1); ++ip )
	{
	    TEST_EQUALITY( ip_coords(cell,ip,0), 1.0*(cell) + 0.5 );
	    TEST_EQUALITY( ip_coords(cell,ip,1), 0.5 );
	    TEST_EQUALITY( ip_coords(cell,ip,2), 0.5 );
	}

	TEST_EQUALITY( dof_val*cell_measure(cell), integrals(cell) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( IntrepidCell, side_ffupdate_all_cells_test )
{
    // INITIALIZATION
    // --------------

    // Basic problem info.
    int dimension = 3;
    int num_elements = 5;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();

    // Build the mesh
    Teuchos::Array<double> element_coordinates =
	buildMeshNodeCoordinates( comm, dimension, num_elements );

    // Build an cell manager for the mesh.
    shards::CellTopology element_topo =
	shards::getCellTopologyData<shards::Hexahedron<8> >();
    shards::CellTopology side_topo =
	shards::getCellTopologyData<shards::Quadrilateral<4> >();
    int degree = 1;
    int side_id = 0;
    DataTransferKit::IntrepidSideCell
	intrepid_cell( side_topo, side_id, element_topo, degree);

    // For each element in the mesh compute cell measures, integrals, and
    // physical cubature points and check them.
    int workset_size = num_elements;
    int num_ip = 1;
    int num_nodes = element_topo.getNodeCount();
    Teuchos::Array<int> coord_dims( 3 );
    coord_dims[0] = workset_size;
    coord_dims[1] = num_nodes;
    coord_dims[2] = dimension;
    Intrepid::FieldContainer<double> element_coords( 
	coord_dims, element_coordinates() );
    Intrepid::FieldContainer<double> cell_measure(workset_size);
    Intrepid::FieldContainer<double> ip_coords;
    Intrepid::FieldContainer<double> integrals(workset_size);
    Intrepid::FieldContainer<double> dofs(workset_size,num_ip);
    int dof_val = 2.0;
    dofs.initialize( dof_val );

    DataTransferKit::IntrepidCell::updateState(
	intrepid_cell, element_coords );

    intrepid_cell.getCellMeasures( cell_measure );
    intrepid_cell.getPhysicalIntegrationCoordinates( ip_coords );
    intrepid_cell.integrate( dofs, integrals );

    TEST_EQUALITY( workset_size, intrepid_cell.getNumCells() );
    TEST_EQUALITY( num_ip, intrepid_cell.getNumIntegrationPoints() );
    TEST_EQUALITY( dimension, intrepid_cell.getSpatialDimension() );

    for ( int cell = 0; cell < num_elements; ++cell )
    {
	for ( int ip = 0; ip < ip_coords.dimension(1); ++ip )
	{
	    TEST_EQUALITY( ip_coords(cell,ip,0), 1.0*(cell) + 0.5 );
	    TEST_EQUALITY( ip_coords(cell,ip,1), 0.0 );
	    TEST_EQUALITY( ip_coords(cell,ip,2), 0.5 );
	}

	TEST_EQUALITY( dof_val*cell_measure(cell), integrals(cell) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( IntrepidCell, point_in_hex_test )
{
    // INITIALIZATION
    // --------------

    // Basic problem info.
    int dimension = 3;
    int num_elements = 5;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();

    // Build the mesh
    Teuchos::Array<double> element_coordinates =
	buildMeshNodeCoordinates( comm, dimension, num_elements );

    // Get the element block part.
    shards::CellTopology element_topo =
	shards::getCellTopologyData<shards::Hexahedron<8> >();

    // Perform point in element checks.
    Intrepid::FieldContainer<double> point(1,dimension);
    point(0,0) = 2.59;
    point(0,1) = 0.54;
    point(0,2) = 0.32;
    DataTransferKit::IntrepidCell hex_cell( element_topo, 1 );
    Teuchos::Array<int> coord_dims( 3 );
    coord_dims[0] = 1;
    coord_dims[1] = 8;
    coord_dims[2] = dimension;
    bool point_in_element = false;
    for ( int elem_id = 0; elem_id < num_elements; ++elem_id )
    {
	Intrepid::FieldContainer<double> element_nodes(
	    coord_dims, element_coordinates(8*dimension*elem_id,dimension*8) );
	hex_cell.allocateCellState( element_nodes );
	point_in_element = hex_cell.pointInPhysicalCell( point, 1.0e-8 );

	if ( 2 == elem_id )
	{
	    TEST_ASSERT( point_in_element );
	}
	else
	{
	    TEST_ASSERT( !point_in_element );
	}
    }
}

//---------------------------------------------------------------------------//
// end tstIntrepidCell
//---------------------------------------------------------------------------//
