//---------------------------------------------------------------------------//
/*!
 * \file MoabMesh.cpp
 * \author Stuart R. Slattery
 * \brief Moab mesh definition more example.
 */
//---------------------------------------------------------------------------//

#include "MoabMesh.hpp"

#include <DTK_Exception.hpp>
#include <DTK_MeshTypes.hpp>

#include <MBCore.hpp>

//---------------------------------------------------------------------------//
MoabMesh::MoabMesh( const RCP_Comm& comm, 
		    const std::string& filename,
		    const moab::EntityType& block_topology,
		    const int partitioning_type )
    : d_comm( comm )
    , d_topology( block_topology )
    , d_node_dim( 0 )
    , d_nodes_per_element( 0 )
{ 
    // Compute the node dimension.
    if ( d_topology == moab::MBTRI )
    {
	d_node_dim = 2;
    }
    else if ( d_topology == moab::MBQUAD )
    {
	d_node_dim = 2;
    }
    else if ( d_topology == moab::MBTET )
    {
	d_node_dim = 3;
    }
    else if ( d_topology == moab::MBHEX )
    {
	d_node_dim = 3;
    }
    else if ( d_topology == moab::MBPYRAMID )
    {
	d_node_dim = 3;
    }
    else
    {
	d_node_dim = 0;
    }

    // Create a moab instance.
    moab::ErrorCode error;
    d_moab = Teuchos::rcp( new moab::Core() );

    // Load the mesh.
    d_moab->load_mesh( &filename[0] );
    moab::EntityHandle root_set = d_moab->get_root_set();

    // Extract the elements with this block's topology.
    std::vector<moab::EntityHandle> global_elements;
    error = d_moab->get_entities_by_type( 
	root_set, d_topology, global_elements );
    DataTransferKit::testInvariant( error == moab::MB_SUCCESS,
				    "Error getting global block elements" );

    // Partition the mesh.
    int my_rank = d_comm->getRank();
    std::vector<moab::EntityHandle> elem_vertices;
    error = d_moab->get_adjacencies( &global_elements[0],
				     1,
				     0,
				     false,
				     elem_vertices );
    DataTransferKit::testInvariant( error == moab::MB_SUCCESS,
				    "Error getting adjacencies." );
    d_nodes_per_element = elem_vertices.size();

    std::vector<double> elem_coords( 3*d_nodes_per_element );
    std::vector<moab::EntityHandle>::const_iterator global_elem_iterator;
    for ( global_elem_iterator = global_elements.begin();
	  global_elem_iterator != global_elements.end();
	  ++global_elem_iterator )
    {
	error = d_moab->get_adjacencies( &*global_elem_iterator,
					 1,
					 0,
					 false,
					 elem_vertices );
	DataTransferKit::testInvariant( error == moab::MB_SUCCESS,
					"Error getting element vertices" );

	error = d_moab->get_coords( &elem_vertices[0], 
				    elem_vertices.size(),
				    &elem_coords[0] );
	DataTransferKit::testInvariant( error == moab::MB_SUCCESS,
					"Error getting element vertices" );

	if ( partitioning_type == 0 )
	{
	    if ( my_rank == 0 )
	    {
		if ( elem_coords[0] <= -2.5 )
		{
		    d_elements.push_back( *global_elem_iterator );
		}
	    }
	    else if ( my_rank == 1 )
	    {
		if ( elem_coords[0] >= -2.5 && elem_coords[0] <= 0.0 )
		{
		    d_elements.push_back( *global_elem_iterator );
		}
	    }
	    else if ( my_rank == 2 )
	    {
		if ( elem_coords[0] >= 0.0 && elem_coords[0] <= 2.5 )
		{
		    d_elements.push_back( *global_elem_iterator );
		}
	    }
	    else
	    {
		if ( elem_coords[0] >= 2.5 )
		{
		    d_elements.push_back( *global_elem_iterator );
		}
	    }
	}
	else if ( partitioning_type == 1 )
	{
	    if ( my_rank == 0 )
	    {
		if ( elem_coords[1] <= -2.5 )
		{
		    d_elements.push_back( *global_elem_iterator );
		}
	    }
	    else if ( my_rank == 1 )
	    {
		if ( elem_coords[1] >= -2.5 && elem_coords[0] <= 0.0 )
		{
		    d_elements.push_back( *global_elem_iterator );
		}
	    }
	    else if ( my_rank == 2 )
	    {
		if ( elem_coords[1] >= 0.0 && elem_coords[0] <= 2.5 )
		{
		    d_elements.push_back( *global_elem_iterator );
		}
	    }
	    else
	    {
		if ( elem_coords[1] >= 2.5 )
		{
		    d_elements.push_back( *global_elem_iterator );
		}
	    }
	}
	else
	{
	    throw DataTransferKit::InvariantException( 
		"Partitioning type not supported." );
	}
    }
    d_comm->barrier();

    // Get the nodes.
    error = d_moab->get_connectivity( &d_elements[0],
				      d_elements.size(),
				      d_vertices );
    DataTransferKit::testInvariant( error == moab::MB_SUCCESS,
				    "Error getting connecting vertices" );

    // Get the node coordinates.
    d_coords.resize( 3 * d_vertices.size() );
    std::vector<double> interleaved_coords( d_coords.size() );
    error = d_moab->get_coords( &d_vertices[0], d_vertices.size(),
				&interleaved_coords[0] );
    DataTransferKit::testInvariant( error == moab::MB_SUCCESS,
				    "Error getting mesh coordinates" );
    for ( int n = 0; n < (int) d_vertices.size(); ++n )
    {
	for ( int d = 0; d < (int) d_node_dim; ++d )
	{
	    d_coords[ d*d_vertices.size() + n ] =
		interleaved_coords[ n*d_node_dim + d ];
	}
    }

    // Get the connectivity.
    int connectivity_size = d_elements.size() * d_nodes_per_element;
    d_connectivity.resize( connectivity_size );
    std::vector<moab::EntityHandle> elem_conn;
    for ( int i = 0; i < (int) d_elements.size(); ++i )
    {
	error = d_moab->get_connectivity( &d_elements[i], 1, elem_conn );

	DataTransferKit::testInvariant( error == moab::MB_SUCCESS,
					"Error getting element vertices." );
	DataTransferKit::testInvariant( elem_conn.size() == d_nodes_per_element,
					"Num adj != nodes per element." );

	for ( int n = 0; n < (int) elem_conn.size(); ++n )
	{
	    d_connectivity[ n*d_elements.size() + i ] = elem_conn[n];
	}
    }

    // Get the permutation vector.
    d_permutation_list.resize( d_nodes_per_element );
    for ( int i = 0; i < (int) d_nodes_per_element; ++i )
    {
	d_permutation_list[i] = i;
    }
}

//---------------------------------------------------------------------------//
void MoabMesh::tag( const ArrayField& data )
{
    moab::ErrorCode error;
    moab::Tag tag;
    error = d_moab->tag_get_handle( "data", 1, moab::MB_TYPE_DOUBLE, tag );
    DataTransferKit::testInvariant( error == moab::MB_SUCCESS,
				    "Error creating data tag." );
    error = d_moab->tag_set_data( tag, 
				  &d_vertices[0], 
				  d_vertices.size(),
				  &data[0] );
    DataTransferKit::testInvariant( error == moab::MB_SUCCESS,
				    "Error tagging data." );
}

//---------------------------------------------------------------------------//
void MoabMesh::write( const std::string& filename )
{
    std::string filename0 = "0" + filename;
    std::string filename1 = "1" + filename;
    std::string filename2 = "2" + filename;
    std::string filename3 = "3" + filename;

    if ( d_comm->getRank() == 0 )
    {
	d_moab->write_mesh( &filename0[0] );
    }
    else if ( d_comm->getRank() == 1 )
    {
	d_moab->write_mesh( &filename1[0] );
    }
    else if ( d_comm->getRank() == 2 )
    {
	d_moab->write_mesh( &filename2[0] );
    }
    else if ( d_comm->getRank() == 3 )
    {
	d_moab->write_mesh( &filename3[0] );
    }
    d_comm->barrier();
}

//---------------------------------------------------------------------------//
std::size_t MoabMesh::blockTopology() const
{ 
    if ( d_topology == moab::MBTRI )
    {
	return DataTransferKit::DTK_TRIANGLE;
    }
    else if ( d_topology == moab::MBQUAD )
    {
	return DataTransferKit::DTK_QUADRILATERAL;
    }
    else if ( d_topology == moab::MBTET )
    {
	return DataTransferKit::DTK_TETRAHEDRON;
    }
    else if ( d_topology == moab::MBHEX )
    {
	return DataTransferKit::DTK_HEXAHEDRON;
    }
    else if ( d_topology == moab::MBPYRAMID )
    {
	return DataTransferKit::DTK_PYRAMID;
    }
    return 0;
}

//---------------------------------------------------------------------------//
// end MoabMesh.cpp
//---------------------------------------------------------------------------//

