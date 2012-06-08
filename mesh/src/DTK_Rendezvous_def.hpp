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
#include <multimap>
#include <algorithm>
#include <cassert>

#include "DTK_Rendezvous.hpp"
#include <DTK_InverseComm.hpp>
#include <DTK_Exception.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Hashtable.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include <Tpetra_Export.hpp>
#include <Tpetra_Vector.hpp>

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
    getMeshInBox( mesh, d_global_box, nodes_in_box, elements_in_box );
        
    // Construct the rendezvous decomposition of the mesh with RCB.
    d_rcb = Teuchos::rcp( new RCB<Mesh>( mesh, nodes_in_box, d_comm ) );
    testPostcondition( d_rcb != Teuchos::null,
		       "Error creating RCB decomposition." );
    d_rcb->partition();

    // Send the mesh to the rendezvous decomposition and build the concrete
    // mesh.
    sendMeshToRendezvous( mesh, elements_in_box );

    // Clear the extracted mesh information.
    nodes_in_box.clear();
    elements_in_box.clear();
    
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
 * \brief Extract the mesh nodes and elements that are in a bounding box.
 */
template<typename Mesh>
void Rendezvous<Mesh>::getMeshInBox( const Mesh& mesh,
				     const BoundingBox& box,
				     std::vector<char>& nodes_in_box,
				     std::vector<char>& elements_in_box )
{
    // Create a map indexed by handle containing the actual handle
    // location. This will give us logarithmic time access to connectivity. I
    // should write a more general hash table to improve this access time.
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
	    nodes_in_box.push_back( box.pointInBox( node_coords ) );
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
	    nodes_in_box.push_back( box.pointInBox( node_coords ) );
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
		node_index = 
		    node_indices.find( *connectivity_iterator )->second;
		nodes_in_box[ node_index ] = 1;
	    }
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send the mesh to the rendezvous decomposition and build the concrete
 * mesh. 
 */
template<typename Mesh>
void Rendezvous<Mesh>sendMeshToRendezvous( 
    const Mesh& mesh, const std::vector<char>& elements_in_box )
{
    // Setup data structures for pulling in mesh data in the rendezvous
    // decomposition.
    std::set<handle_type> rendezvous_nodes;
    std::vector<double> rendezvous_coords;
    int element_type;
    int element_topology;
    int nodes_per_element;
    std::set<handle_type> rendezvous_elements;
    std::vector<handle_type> rendezvous_connectivity;

    // Setup the communication patterns for moving the mesh to the rendezvous
    // decomposition. This will also move the node and element handles to the
    // rendezvous decomposition.
    RCP_TpetraMap export_node_map, import_node_map;
    RCP_TpetraMap export_element_map, import_element_map;
    setupCommunication( mesh, elements_in_box,
			export_node_map, import_node_map,
			export_element_map, import_element_map );
    Tpetra::Export<handle_type> node_exporter( export_node_map, 
					       import_node_map );
    Tpetra::Export<handle_type> element_exporter( export_element_map, 
						  import_element_map );

    // Get the rendezvous nodes and elements out of the maps.

    // Move the node coordinates to the rendezvous decomposition. 
    Teuchos::RCP< Tpetra::vector<double> > export_x_coords;
    Teuchos::RCP< Tpetra::vector<double> > export_y_coords;
    Teuchos::RCP< Tpetra::vector<double> > export_z_coords;

    Teuchos::ArrayRCP<double> x_coords_view;
    Teuchos::ArrayRCP<double> y_coords_view;
    Teuchos::ArrayRCP<double> z_coords_view;

    MT::const_coordinate_iterator export_coords = MT::coordsBegin( mesh );
    if ( MT::interleavedCoordinates( mesh ) )
    {
	x_coords_view = Teuchos::ArrayRCP<double>( num_nodes );
	y_coords_view = Teuchos::ArrayRCP<double>( num_nodes );
	z_coords_view = Teuchos::ArrayRCP<double>( num_nodes );
	Teuchos::ArrayRCP<double>::iterator x_iterator;
	Teuchos::ArrayRCP<double>::iterator y_iterator;
	Teuchos::ArrayRCP<double>::iterator z_iterator;
	for ( x_iterator = x_coords_view.begin(),
	      y_iterator = y_coords_view.begin(),
	      z_iterator = z_coords_view.begin();
	      x_iterator != x_coords_view.end();
	      ++x_iterator, ++y_iterator, ++z_iterator )
	{
	    *x_iterator = *export_coords;
	    ++export_coords;
	    *y_iterator = *export_coords;
	    ++export_coords;
	    *z_iterator = *export_coords;
	    ++export_coords;
	}
    }
    else
    {
	x_coords_view = Teuchos::ArrayRCP<double>( 
	    &*export_coords, 0, num_nodes, false );
	y_coords_view = Teuchos::ArrayRCP<double>( 
	    &(*export_coords + num_nodes), 0, num_nodes, false );
	z_coords_view = Teuchos::ArrayRCP<double>( 
	    &(*export_coords + 2*num_nodes), 0, num_nodes, false );
    }

    Teuchos::RCP< Tpetra::vector<double> > import_x_coords = 
	createVector( import_node_map );
    Teuchos::RCP< Tpetra::vector<double> > import_y_coords =
	createVector( import_node_map );
    Teuchos::RCP< Tpetra::vector<double> > import_z_coords =
	createVector( import_node_map );

    import_x_coords->doExport( 
	*export_x_coords, node_exporter, Tpetra::Insert );
    import_y_coords->doExport( 
	*export_y_coords, node_exporter, Tpetra::Insert );
    import_z_coords->doExport( 
	*export_z_coords, node_exporter, Tpetra::Insert );

    // Free up some memory now that we've moved the node data over.
    x_coords_view.clear();
    y_coords_view.clear();
    z_coords_view.clear();

    // Move the element connectivity to the rendezvous decomposition.

    // Construct the mesh container from the collected data, effectively
    // wrapping it with mesh traits. This may not be a good idea. I may create
    // another mesh creation function to directly accept these Tpetra vectors
    // instead as this is making a temporary copy.
    MeshContainer mesh_container( rendezvous_nodes, 
				  rendezvous_coords,
				  MT::elementType( mesh ),
				  MT::elementTopology( mesh ),
				  MT::nodesPerElement( mesh ),
				  rendezvous_elements, 
				  rendezvous_connectivity );

    // Clear the collected data.
    rendezvous_nodes.clear();
    rendezvous_coords.clear();
    rendezvous_elements.clear();
    rendezvous_connectivity.clear();

    // Build the rendezvous mesh from the mesh container.
    d_rendezvous_mesh = createRendezvousMesh( mesh_container );

    testPostcondition( d_rendezvous_mesh != Teuchos::null,
		       "Error creating rendezvous mesh." );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Setup the communication patterns.
 */
template<typename Mesh>
void Rendezvous<Mesh>::setupCommunication( 
    const Mesh& mesh,
    const std::vector<char>& elements_in_box,
    RCP_TpetraMap& export_node_map, 
    RCP_TpetraMap& import_node_map,
    RCP_TpetraMap& export_element_map, 
    RCP_TpetraMap& import_element_map )
{
    // Start by building the Tpetra map in the original decomposition for the
    // nodes.
    int num_nodes = std::distance( MT::nodesBegin( mesh ), 
				   MT::nodesEnd( mesh ) );
    Teuchos::ArrayView<handle_type> export_node_view( 
	&*MT::nodesBegin( mesh ), num_nodes );
    Teuchos::RCP< Tpetra::Map<handle_type> > export_node_map = 
	Tpetra::createNonContigMap<handle_type>( export_node_view, 
						 d_global_comm );
    testPostcondition( export_node_map != Teuchos::null,
		       "Error creating node export map." );

    // Do the same thing for the elements in the original decomposition.
    int num_elements = std::distance( MT::elementsBegin( mesh ), 
				   MT::elementsEnd( mesh ) );
    Teuchos::ArrayView<handle_type> export_element_view( 
	&*MT::elementsBegin( mesh ), num_elements );
    Teuchos::RCP< Tpetra::Map<handle_type> > export_element_map = 
	Tpetra::createNonContigMap<handle_type>( export_element_view, 
						 d_global_comm );
    testPostcondition( export_element_map != Teuchos::null,
		       "Error creating element export map." );

    // Create an index map for logarithmic time access to connectivity data.
    std::map<handle_type,int> node_indices;
    for ( export_node_iterator = MT::nodesBegin( mesh ),
		  int n = 0;
	  export_node_iterator != MT::nodesEnd( mesh );
	  ++export_node_iterator, ++n )
    {
	node_indices[ *node_iterator ] = n;
    }

    // Create a multimap of the destination procs for all local elements and
    // connecting nodes in the global bounding box.
    std::multimap<int,handle_type> node_proc_map;
    std::multimap<int,handle_type> element_proc_map;
    int nodes_per_element = MT::nodesPerElement( mesh );
    int node_index;
    int destination_proc;
    double node_coords[3];
    MT::const_handle_iterator mesh_nodes = MT::nodesBegin( mesh );
    MT::const_coordinate_iterator mesh_coords = MT::coordsBegin( mesh );
    MT::const_handle_iterator element_iterator;
    MT::const_handle_iterator connectivity_iterator;
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
		node_index = 
		    node_indices.find( *connectivity_iterator )->second;
		node_coords[0] = mesh_coords[ 3*node_index ];
		node_coords[1] = mesh_coords[ 3*node_index + 1 ];
		node_coords[2] = mesh_coords[ 3*node_index + 2 ];
		destination_proc = d_rcb->getDestinationProc( node_coords );
		node_proc_map.insert( 
		    std::pair<int,handle_type>(
			destination_proc, mesh_nodes[ node_index ] ) );
		element_proc_map.insert( 
		    std::pair<int,handle_type>( 
			destination_proc, *element_iterator ) );
	    }
	}
    }

    // Now we know where the nodes and elements need to go. Build the maps
    // through an inverse communication operation.
    InverseComm inverse_comm( d_global_comm );
    import_node_map = inverse_comm.createImportMap( node_proc_map );
    import_element_map = inverse_comm.createImportMap( element_proc_map );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_RENDEZVOUS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_Rendezvous_def.hpp
//---------------------------------------------------------------------------//
