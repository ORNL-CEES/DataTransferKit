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

#include "DTK_Rendezvous.hpp"
#include <DTK_Exception.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Tpetra_Distributor.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Mesh>
Rendezvous<Mesh>::Rendezvous( const RCP_Comm& global_comm,
			      const BoundingBox& global_box )
    : d_global_comm( global_comm )
    , d_global_box( global_box )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Mesh>
Rendezvous<Mesh>::~Rendezvous()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Build the rendezvous decomposition.
 */
template<class Mesh> 
void Rendezvous<Mesh>::build( const Mesh& mesh )
{
    // Extract the mesh nodes and elements that are in the bounding box.
    std::vector<char> nodes_in_box;
    std::vector<char> elements_in_box;
    getMeshInBox( mesh, d_global_box, nodes_in_box, elements_in_box );
        
    // Construct the rendezvous decomposition of the mesh with RCB using the
    // nodes that are in the box.
    d_rcb = Teuchos::rcp( new RCB<Mesh>( mesh, nodes_in_box, d_comm ) );
    testPostcondition( d_rcb != Teuchos::null,
		       "Error creating RCB decomposition." );
    d_rcb->partition();

    // Send the mesh in the box to the rendezvous decomposition and build the
    // concrete mesh.
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
template<class Mesh>
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
template<class Mesh>
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
template<class Mesh>
void Rendezvous<Mesh>::getMeshInBox( const Mesh& mesh,
				     const BoundingBox& box,
				     std::vector<char>& nodes_in_box,
				     std::vector<char>& elements_in_box )
{
    // Create a map indexed by node handle containing the actual node handle
    // location. This will give us logarithmic time access to connectivity. I
    // should write a more general hash table to improve this access time as
    // I'm using this strategy for most mesh operations.
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

    // Get all of the nodes that are in the box. 
    double node_coords[3];
    typename MT::const_coordinate_iterator coord_iterator;
    for ( coord_iterator = MT::coordsBegin( mesh );
	  coord_iterator != MT::coordsEnd( mesh ); )
    {
	for ( int i = 0; i < 3; ++i, ++coord_iterator )
	{
	    node_coords[i] = *coord_iterator;
	}
	nodes_in_box.push_back( box.pointInBox( node_coords ) );
    }
    assert( nodes_in_box.size() == num_nodes );

    // For those nodes that are in the box, get the elements that they
    // construct. These elements are in the box.
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
	this_element_in_box = 0;
	for ( int n = 0; n < nodes_per_element; ++n, ++connectivity_iterator )
	{
	    node_index = node_indices.find( *connectivity_iterator )->second;
	    if ( nodes_in_box[ node_index ] )
	    {
		this_element_in_box = 1;
	    }
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
	if ( elements_in_box[i] )
	{
	    for ( int n = 0; n < nodes_per_element; ++n, ++connectivity_iterator )
	    {
		node_index = 
		    node_indices.find( *connectivity_iterator )->second;
		nodes_in_box[ node_index ] = 1;
	    }
	}
	else
	{
	    connectivity_iterator += nodes_per_element;
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send the mesh to the rendezvous decomposition and build the concrete
 * mesh. 
 */
template<class Mesh>
void Rendezvous<Mesh>sendMeshToRendezvous( 
    const Mesh& mesh, const std::vector<char>& elements_in_box )
{
    // Setup the communication patterns for moving the mesh to the rendezvous
    // decomposition. This will also move the node and element handles to the
    // rendezvous decomposition.
    std::set<handle_type> rendezvous_nodes;
    std::set<handle_type> rendezvous_elements;
    setupCommunication( mesh, elements_in_box,
			rendezvous_nodes, rendezvous_elements );

    // Setup export node map.
    int num_nodes = std::distance( MT::nodesBegin( mesh ), 
				   MT::nodesEnd( mesh ) );
    Teuchos::ArrayView<handle_type> export_node_view( 
	&*MT::nodesBegin( mesh ), num_nodes );
    RCP_TpetraMap export_node_map = Tpetra::createNonContigMap<handle_type>( 
	export_node_view, d_global_comm );
    testPostcondition( export_node_map != Teuchos::null,
		       "Error creating node export map." );

    // Setup import node map.
    Teuchos::ArrayView<handle_type> rendezvous_nodes_view( 
	&*rendezvous_nodes.begin(), rendezvous_nodes.size() );
    RCP_TpetraMap import_node_map = Tpetra::createNonContigMap<handle_type>(
	import_node_view, d_global_comm );
    testPostcondition( import_node_map != Teuchos::null,
		       "Error creating node import map." );

    // Setup export element map.
    int num_elements = std::distance( MT::elementsBegin( mesh ), 
				      MT::elementsEnd( mesh ) );
    Teuchos::ArrayView<handle_type> export_element_view( 
	&*MT::elementsBegin( mesh ), num_elements );
    RCP_TpetraMap export_element_map = Tpetra::createNonContigMap<handle_type>(
	export_element_view, d_global_comm );
    testPostcondition( export_element_map != Teuchos::null,
		       "Error creating element export map." );
   
    // Setup import element map.
    Teuchos::ArrayView<handle_type> rendezvous_elements_view( 
	&*rendezvous_elements.begin(), rendezvous_elements.size() );
    RCP_TpetraMap import_element_map = Tpetra::createNonContigMap<handle_type>(
	import_element_view, d_global_comm );
    testPostcondition( import_element_map != Teuchos::null,
		       "Error creating element import map." );

    // Setup exporters.
    Tpetra::Export<handle_type> node_exporter( export_node_map, 
					       import_node_map );
    Tpetra::Export<handle_type> element_exporter( export_element_map, 
						  import_element_map );

    // Move the node coordinates to the rendezvous decomposition.
    int num_coords = 3*num_nodes;
    Teuchos::ArrayRCP<double> export_coords_view( 
	&*MT::coordsBegin( mesh ), num_coords, false );
    Teuchos::RCP< Tpetra::MultiVector<double> > export_coords = 
	createMultiVectorFromView( export_node_map, export_coords_view, 1, 3 );
    Tpetra::MultiVector<double> import_coords( import_node_map, 3 );
    import_coords.doExport( *export_coords, node_exporter, Tpetra::Insert );

    // Move the element connectivity to the rendezvous decomposition. 
    int num_conn = MT::nodesPerElement( mesh ) * num_elements;
    Teuchos::ArrayRCP<handle_type> export_conn_view( 
	&*MT::connectivityBegin( mesh ), num_conn, false );
    Teuchos::RCP< Tpetra::MultiVector<handle_type> > export_conn = 
	createMultiVectorFromView( export_element_map, export_conn_view, 
				   1, MT::nodesPerElement( mesh ) );
    Tpetra::MultiVector<handle_type> import_conn( import_node_map, 
						  MT::nodesPerElement( mesh ) );
    import_conn.doExport( *export_conn, node_exporter, Tpetra::Insert );

    // Construct the mesh container from the collected data, effectively
    // wrapping it with mesh traits. 
    MeshContainer mesh_container( rendezvous_nodes, 
				  import_coords.get1dView(),
				  MT::elementType( mesh ),
				  MT::elementTopology( mesh ),
				  MT::nodesPerElement( mesh ),
				  rendezvous_elements, 
				  import_conn.get1dView() );

    // Clear the collected data.
    rendezvous_nodes.clear();
    rendezvous_elements.clear();

    // Build the concrete rendezvous mesh from the mesh container.
    d_rendezvous_mesh = createRendezvousMesh( mesh_container );
    testPostcondition( d_rendezvous_mesh != Teuchos::null,
		       "Error creating rendezvous mesh." );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Setup the communication patterns.
 */
template<class Mesh>
void Rendezvous<Mesh>::setupCommunication( 
    const Mesh& mesh,
    const std::vector<char>& elements_in_box,
    std::set<handle_type>& rendezvous_nodes,
    std::set<handle_type>& rendezvous_elements )
{
    // Create a node index map for logarithmic time access to connectivity
    // data. 
    std::map<handle_type,int> node_indices;
    for ( export_node_iterator = MT::nodesBegin( mesh ),
			 int n = 0;
	  export_node_iterator != MT::nodesEnd( mesh );
	  ++export_node_iterator, ++n )
    {
	node_indices[ *node_iterator ] = n;
    }

    // Create a element index map for logarithmic time access to connectivity
    // data. 
    std::map<handle_type,int> element_indices;
    for ( export_element_iterator = MT::elementsBegin( mesh ),
			    int n = 0;
	  export_element_iterator != MT::elementsEnd( mesh );
	  ++export_element_iterator, ++n )
    {
	element_indices[ *element_iterator ] = n;
    }

    // Get destination procs for all local elements in the global bounding
    // box. The element will need to be sent to each partition that its
    // connecting nodes exist in. We'll make a unique destination proc set for
    // each element.
    std::vector< std::set<int> > export_element_procs_set( elements_in_box.size() );
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
	if ( elements_in_box[i] )
	{
	    for ( int n = 0; n < nodes_per_element; ++n, ++connectivity_iterator )
	    {
		node_index = 
		    node_indices.find( *connectivity_iterator )->second;
		node_coords[0] = mesh_coords[ 3*node_index ];
		node_coords[1] = mesh_coords[ 3*node_index + 1 ];
		node_coords[2] = mesh_coords[ 3*node_index + 2 ];
		destination_proc = d_rcb->getDestinationProc( node_coords );
		export_element_procs_set[i].insert( destination_proc );
	    }
	}
	else
	{
	    connectivity_iterator += nodes_per_element;
	}
    }

    // Unroll the vector of sets into two vectors; one containing the element
    // handle and the other containing the corresponding element destination.
    std::vector<handle_type> export_elements;
    std::vector<int> export_element_procs;
    std::vector< std::set<int> >::const_iterator element_vec_iterator;
    std::set<int>::const_iterator proc_set_iterator;
    for ( element_vec_iterator = export_elements_procs_set.begin(),
	      element_iterator = MT::beginElements( mesh );
	  element_vec_iterator != export_elements_procs_set.end();
	  ++element_vec_iterator, ++element_iterator )
    {
	for ( proc_set_iterator = element_vec_iterator->begin();
	      proc_set_iterator != element_vec_iterator->end();
	      ++proc_set_iterator )
	{
	    export_elements.push_back( *element_iterator );
	    export_element_procs.push_back( *proc_set_iterator );
	}
    }
    export_elements_procs_set.clear();

    // Now we know where the elements need to go. Move the elements to the
    // rendezvous decomposition through an inverse communciation operation.
    Tpetra::Distributor element_distributor( d_global_comm );
    int num_import_elements = element_distributor.createFromSends(
	Teuchos::arrayViewFromVector( export_element_procs ) );
    std::vector<handle_type> import_elements( num_import_elements );
    element_distributor.doPostsAndWaits( 
	Teuchos::arrayViewFromVector( export_elements ), 1,
	Teuchos::arrayViewFromVector( import_elements ) );

    // Next move these into the rendezvous element set so that we have a
    // unique list of the elements.
    std::vector<handle_type>::const_iterator import_element_iterator;
    for ( import_element_iterator = import_elements.begin();
	  import_element_iterator != import_elements.end();
	  ++import_element_iterator )
    {
	rendezvous_elements.insert( *import_element_iterator );
    }
    import_elements.clear();

    // Now get the destination procs for all the nodes. This will be the same
    // destination procs as all of their parent elements. Therefore, nodes may
    // then also have to go to multiple procs because of this and these procs
    // may be different than their original RCB procs.
    int node_handle;
    MT::const_connectivity_iterator export_conn = MT::connectivityBegin( mesh );
    std::vector<handle_type>::const_iterator export_elements_iterator;
    std::vector<int>::const_iterator export_element_procs_iterator;
    std::vector< std::set<int> > export_node_procs_set;
    for ( export_elements_iterator = export_elements.begin(),
     export_element_procs_iterator = export_element_procs.begin();
	  export_elements_iterator != export_elements.end();;
	  ++export_elements_iterator, ++export_element_procs_iterator )
    {
	element_index = 
	    element_indices.find( *export_elements_iterator )->second;

	for ( int n = 0; n < nodes_per_element; ++n )
	{
	    node_handle = export_conn[ element_index*nodes_per_element + n ];

	    node_index = 
		node_indices.find( node_handle )->second;
	    
	    export_nodes_procs_set[ node_index ].insert( 
		*export_element_procs_iterator );
	}
    }
    export_elements.clear();    
    export_element_procs.clear();
    node_indices.clear();
    element_indices.clear();

    // Unroll the vector of sets into two vectors; one containing the node
    // handle and the other containing the corresponding node destination.
    std::vector<handle_type> export_nodes;
    std::vector<int> export_node_procs;
    std::vector< std::set<int> >::const_iterator node_vec_iterator;
    std::set<int>::const_iterator proc_set_iterator;
    for ( node_vec_iterator = export_nodes_procs_set.begin(),
	      node_iterator = MT::beginNodes( mesh );
	  node_vec_iterator != export_nodes_procs_set.end();
	  ++node_vec_iterator, ++node_iterator )
    {
	for ( proc_set_iterator = node_vec_iterator->begin();
	      proc_set_iterator != node_vec_iterator->end();
	      ++proc_set_iterator )
	{
	    export_nodes.push_back( *node_iterator );
	    export_node_procs.push_back( *proc_set_iterator );
	}
    }
    export_nodes_procs_set.clear();

    // Now we know where the nodes need to go. Move the nodes to the
    // rendezvous decomposition through an inverse communciation operation.
    Tpetra::Distributor node_distributor( d_global_comm );
    int num_import_nodes = node_distributor.createFromSends(
	Teuchos::arrayViewFromVector( export_node_procs ) );
    std::vector<handle_type> import_nodes( num_import_nodes );
    node_distributor.doPostsAndWaits( 
	Teuchos::arrayViewFromVector( export_nodes ), 1,
	Teuchos::arrayViewFromVector( import_nodes ) );
    export_nodes.clear();
    export_node_procs.clear();

    // Next move these into the rendezvous node set so that we have a unique
    // list of the nodes.
    std::vector<handle_type>::const_iterator import_node_iterator;
    for ( import_node_iterator = import_nodes.begin();
	  import_node_iterator != import_nodes.end();
	  ++import_node_iterator )
    {
	rendezvous_nodes.insert( *import_node_iterator );
    }
    import_nodes.clear();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_RENDEZVOUS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_Rendezvous_def.hpp
//---------------------------------------------------------------------------//
