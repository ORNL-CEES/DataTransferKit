//---------------------------------------------------------------------------//
/*!
 * \file DTK_Rendezvous_def.hpp
 * \author Stuart R. Slattery
 * \brief Rendezvous definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_RENDEZVOUS_DEF_HPP
#define DTK_RENDEZVOUS_DEF_HPP

#include <map>
#include <algorithm>
#include <cassert>

#include <DTK_Exception.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ENull.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename Mesh>
Rendezvous<Mesh>::Rendezvous( const RCP_Comm& global_comm,
			      const BoundingBox& global_box )
    : d_global_comm( global_comm )
    , d_global_box( global_box )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<typename Mesh>
Rendezvous<Mesh>::~Rendezvous()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Build the rendezvous decomposition.
 */
template<typename Mesh> 
void Rendezvous<Mesh>::build( const Mesh& mesh )
{
    // Extract the mesh nodes and elements that are in the bounding box.
    std::vector<char> nodes_in_box;
    std::vector<char> elements_in_box;
    getMeshInBox( mesh, nodes_in_box, elements_in_box );
        
    // Construct the rendezvous decomposition of the mesh with RCB.
    d_rcb = Teuchos::rcp( new RCB<Mesh>( mesh, nodes_in_box, d_comm ) );
    testPostcondition( d_rcb != Teuchos::null,
		       "Error creating RCB decomposition." );
    d_rcb->partition();

    // Send the mesh to the rendezvous decomposition.
    sendMeshToRendezvous();

    // Clear the extracted mesh information.
    nodes_in_box.clear();
    elements_in_box.clear();

    // Build the concrete mesh database in the rendezvous decomposition.
    d_rendezvous_mesh = createRendezvousMesh( mesh_container );
    testPostcondition( d_rendezvous_mesh != Teuchos::null,
		       "Error creating rendezvous mesh." );
    
    // Create a kD-tree in the rendezvous decomposition.
    d_kdtree = Teuchos::rcp( new KDTree<handle_type>( d_rendezvous_mesh ) );
    testPostcondition( d_kdtree != Teuchos::null,
		       "Error creating rendezvous kD-tree." );
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Get the rendezvous processes for a list of points.
 */
template<typename Mesh>
std::vector<int> Rendezvous<Mesh>::getRendezvousProcs(
    const std::vector<double>& coords ) const
{
    testPrecondition( coords.size() % 3 == 0, 
		      "Three dimensional coordinates not provided." );
    
    int num_points = coords.size() / 3;
    std::vector<int> destination_procs;
    for ( int i = 0; i < num_points; ++i )
    {
	double point[3] = { coords[3*i], coords[3*i+1], coords[3*i+2] };
	int rendezvous_proc = d_rcb->getDestinationProc( point );
	destination_procs.push_back( rendezvous_proc );
    }

    testPostcondition( (int) destination_procs.size() == num_points,
		       "Error getting destination processes." );

    return destination_procs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the native mesh elements containing a list of coordinates.
 */
template<typename Mesh>
std::vector< Rendezvous<Mesh>::handle_type >
Rendezvous<Mesh>::getElements( const std::vector<double>& coords ) const
{
    testPrecondition( coords.size() % 3 == 0, 
		      "Three dimensional coordinates not provided." );
    
    int num_points = coords.size() / 3;
    std::vector<handle_type> element_handles;
    for ( int i = 0; i < num_points; ++i )
    {
	double point[3] = { coords[3*i], coords[3*i+1], coords[3*i+2] };
	handle_type handle = d_kdtree->findPoint( point );
	element_handles.push_back( handle );
    }

    testPostcondition( (int) element_handles.size() == num_points,
		       "Error getting mesh elements." );

    return element_handles;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Extract the mesh nodes and elements that are in the bounding box.
 */
template<typename Mesh>
void Rendezvous<Mesh>::getMeshInBox( const Mesh& mesh,
				     std::vector<char>& nodes_in_box,
				     std::vector<char>& elements_in_box )
{
    // Create a map indexed by handle containing the actual handle
    // location. This will give us logarithmic time access to connectivity.
    int num_nodes = std::distance( MT::nodesBegin( mesh ),
				   MT::nodesEnd( mesh ) );
    std::map<handle_type,int> node_indices;
    typename MT::const_handle_iterator node_iterator;
    for ( node_iterator = MT::nodesBegin( mesh ),
		  int n = 0;
	  node_iterator != MT::nodesEnd( mesh );
	  ++node_iterator, ++n )
    {
	node_indices[ *node_iterator ] = n;
    }

    // Get all of the nodes that are in the box. We handle the interleaved
    // case and blocked case slightly differently.
    double node_coords[3];
    typename MT::const_coordinate_iterator coord_iterator;
    if ( MT::interleavedCoordinates( mesh ) )
    {
	for ( coord_iterator = MT::coordsBegin( mesh );
	      coord_iterator != MT::coordsEnd( mesh ); )
	{
	    for ( int i = 0; i < 3; ++i, ++coord_iterator )
	    {
		node_coords[i] = *coord_iterator;
	    }
	    nodes_in_box.push_back( 
		(char) d_global_box.pointInBox( node_coords ) );
	}
    }    
    else
    {
	std::vector<double> interleaved_coords( 3*num_nodes );
	int node, dim;
	for ( coord_iterator = MT::coordsBegin( mesh ),
		       int i = 0;
	      coord_iterator != MT::coordsEnd( mesh );
	      ++coord_iterator, ++i )
	{
	    dim = std::floor( i / num_nodes );
	    node = i - dim*num_nodes;
	    interleaved_coords[ 3*node + dim ] = *coord_iterator;
	}
	std::vector<double>::const_iterator interleaved_iterator;
	for ( interleaved_iterator = interleaved_coords.begin();
	      interleaved_coords != interleaved_coords.end(); )
	{
	    for ( int i = 0; i < 3; ++i, ++interleaved_iterator )
	    {
		node_coords[i] = *interleaved_iterator;
	    }
	    nodes_in_box.push_back( 
		(char) d_global_box.pointInBox( node_coords ) );
	}
    }
    assert( nodes_in_box.size() == num_nodes );

    // For those nodes that are in the box, get the elements that they
    // construct.
    int num_elements = std::distance( MT::elementsBegin( mesh ),
				      MT::elementsEnd( mesh ) );
    int nodes_per_element = MT::nodesPerElement( mesh );
    int node_index;
    char this_element_in_box;
    MT::const_handle_iterator element_iterator;
    MT::const_handle_iterator connectivity_iterator;
    for ( element_iterator = MT::beginElements( mesh ),
     connectivity_iterator = MT::beginConnectivity( mesh );
	  element_iterator != MT::endElements( mesh );
	  ++element_iterator )
    {
	for ( int n = 0; n < nodes_per_element; ++n, ++connectivity_iterator )
	{
	    node_index = node_indices.find( *connectivity_iterator )->second;
	    this_element_in_box = nodes_in_box[ node_index ];
	}
	elements_in_box.push_back( this_element_in_box );
    }
    assert( elements_in_box.size() == num_elements );

    // Get the nodes that belong to the elements in the box, but are not in
    // the box themselves. These will also be used in RCB.
    for ( element_iterator = MT::beginElements( mesh ),
     connectivity_iterator = MT::beginConnectivity( mesh ),
		     int i = 0;
	  element_iterator != MT::endElements( mesh );
	  ++element_iterator, ++i )
    {
	for ( int n = 0; n < nodes_per_element; ++n, ++connectivity_iterator )
	{
	    if ( elements_in_box[i] )
	    {
		node_index = node_indices.find( *connectivity_iterator )->second;
		nodes_in_box[ node_index ] = (char) 1;
	    }
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_RENDEZVOUS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_Rendezvous_def.hpp
//---------------------------------------------------------------------------//

