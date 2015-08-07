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
 * \file   tstIntrepidCellLocalMap.cpp
 * \author Stuart Slattery
 * \brief  IntrepidCellLocalMap tests.
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

#include "DTK_IntrepidCellLocalMap.hpp"

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
    coord[1] = i_plane == 1 || i_plane == 2 ? 2.0 : 0.0 ;
    coord[2] = i_plane == 2 || i_plane == 3 ? 2.0 : 0.0 ;
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
TEUCHOS_UNIT_TEST( IntrepidCellLocalMap, element_single_cell_test )
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

    // Create a cell topology.
    shards::CellTopology element_topo =
	shards::getCellTopologyData<shards::Hexahedron<8> >();

    // For each element in the mesh compute cell measures, and
    // physical cubature points and check them.
    int workset_size = 1;
    int num_nodes = element_topo.getNodeCount();
    Teuchos::Array<int> coord_dims( 3 );
    coord_dims[0] = workset_size;
    coord_dims[1] = num_nodes;
    coord_dims[2] = dimension;
    Teuchos::Array<double> coords( dimension );
    Teuchos::Array<double> param_coords( dimension );
    Teuchos::Array<double> physical_coords( dimension );
    Teuchos::Array<double> centroid( dimension );
    
    for ( int cell = 0; cell < num_elements; ++cell )
    {
	// Create a point.
	coords[0] = 1.0*(cell) + 0.5;
	coords[1] = 0.5;
	coords[2] = 1.5;

	// Create an element.
	Intrepid::FieldContainer<double> element_coords( 
	    coord_dims, element_coordinates(num_nodes*cell*dimension,dimension*num_nodes) );

	// Test the measure.
	double measure = 
	    DataTransferKit::IntrepidCellLocalMap::measure( element_topo, element_coords );
	TEST_EQUALITY( measure, 4.0 );

	// Test the centroid.
	DataTransferKit::IntrepidCellLocalMap::centroid( element_topo,
							 element_coords,
							 centroid );
	TEST_EQUALITY( centroid[0], 1.0*(cell) + 0.5 );
	TEST_EQUALITY( centroid[1], 1.0 );
	TEST_EQUALITY( centroid[2], 1.0 );

	// Test the reference frame map.
	DataTransferKit::IntrepidCellLocalMap::mapToReferenceFrame( element_topo,
								    element_coords,
								    coords,
								    param_coords );
	TEST_EQUALITY( param_coords[0], 0.0 );
	TEST_EQUALITY( param_coords[1], -0.5 );
	TEST_EQUALITY( param_coords[2], 0.5 );

	// Test the point inclusion.
	bool point_inclusion = 
	    DataTransferKit::IntrepidCellLocalMap::checkPointInclusion( element_topo,
								       param_coords,
								       1.0e-6 );
	TEST_ASSERT( point_inclusion );

	// Test the physical frame map.
	DataTransferKit::IntrepidCellLocalMap::mapToPhysicalFrame( element_topo,
								   element_coords,
								   param_coords,
								   physical_coords );
	TEST_EQUALITY( physical_coords[0], coords[0] );
	TEST_EQUALITY( physical_coords[1], coords[1] );
	TEST_EQUALITY( physical_coords[2], coords[2] );
    }
}

//---------------------------------------------------------------------------//
// end tstIntrepidCellLocalMap.cpp
//---------------------------------------------------------------------------//
