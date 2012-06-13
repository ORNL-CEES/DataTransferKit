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

#include "DTK_MeshContainer.hpp"
#include <DTK_Exception.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Tpetra_Distributor.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_MultiVector.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Mesh>
Rendezvous<Mesh>::Rendezvous( const RCP_Comm& comm,
			      const BoundingBox& global_box )
    : d_comm( comm )
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
    Teuchos::Array<int> nodes_in_box;
    Teuchos::Array<int> elements_in_box;
    getMeshInBox( mesh, d_global_box, nodes_in_box, elements_in_box );

    // Construct the rendezvous decomposition of the mesh with RCB using the
    // nodes that are in the box.
    d_rcb = Teuchos::rcp( new RCB<Mesh>( mesh, 
					 Teuchos::arcpFromArray( nodes_in_box ), 
					 d_comm ) );
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
 * \brief Get the rendezvous processes for a list of points. ( Need to specify
 * blocked or interleaved here. Should probably do blocked for constistency.)
 */
template<class Mesh>
Teuchos::Array<int> Rendezvous<Mesh>::getRendezvousProcs(
    const Teuchos::Array<double>& coords ) const
{
    testPrecondition( coords.size() % 3 == 0, 
		      "Three dimensional coordinates not provided." );
    
    double point[3];
    ordinal_type num_points = coords.size() / 3;
    int rendezvous_proc;
    Teuchos::Array<int> destination_procs;
    for ( ordinal_type i = 0; i < num_points; ++i )
    {
	point[0] = coords[ i ];
	point[1] = coords[ num_points + i ];
	point[2] = coords[ 2*num_points + i ];
	rendezvous_proc = d_rcb->getDestinationProc( point );
	destination_procs.push_back( rendezvous_proc );
    }

    testPostcondition( static_cast<ordinal_type>( destination_procs.size() )
		       == num_points,
		       "Error getting destination processes." );

    return destination_procs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the native mesh elements containing a list of coordinates.  (
 * Need to specify blocked or interleaved here. Should probably do blocked for
 * constistency.)
 */
template<class Mesh>
Teuchos::Array< typename Rendezvous<Mesh>::handle_type >
Rendezvous<Mesh>::getElements( const Teuchos::Array<double>& coords ) const
{
    testPrecondition( coords.size() % 3 == 0, 
		      "Three dimensional coordinates not provided." );

    double point[3];
    handle_type element_handle;
    ordinal_type num_points = coords.size() / 3;
    Teuchos::Array<handle_type> element_handles;
    for ( ordinal_type i = 0; i < num_points; ++i )
    {
	point[0] = coords[ i ];
	point[1] = coords[ num_points + i ];
	point[2] = coords[ 2*num_points + i ];
	element_handle = d_kdtree->findPoint( point );
	element_handles.push_back( element_handle );
    }

    testPostcondition( static_cast<ordinal_type>( element_handles.size() )
		       == num_points,
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
				     Teuchos::Array<int>& nodes_in_box,
				     Teuchos::Array<int>& elements_in_box )
{
    // Create a map indexed by node handle containing the actual node handle
    // location. This will give us logarithmic time access to connectivity. I
    // should write a more general hash table to improve this access time as
    // I'm using this strategy for most mesh operations.
    ordinal_type num_nodes = std::distance( MT::nodesBegin( mesh ),
					    MT::nodesEnd( mesh ) );
    std::map<handle_type,ordinal_type> node_indices;
    typename MT::const_node_iterator node_iterator;
    ordinal_type m = 0;
    for ( node_iterator = MT::nodesBegin( mesh );
	  node_iterator != MT::nodesEnd( mesh );
	  ++node_iterator )
    {
	node_indices[ *node_iterator ] = m;
	++m;
    }

    // Get all of the nodes that are in the box. 
    double node_coords[3];
    typename MT::const_coordinate_iterator mesh_coords = 
	MT::coordsBegin( mesh );
    for ( ordinal_type n = 0; n < num_nodes; ++n )
    {
	for ( int i = 0; i < 3; ++i )
	{
	    node_coords[i] = mesh_coords[ i*num_nodes + n ];
	}
	nodes_in_box.push_back( box.pointInBox( node_coords ) );
    }
    assert( (ordinal_type) nodes_in_box.size() == num_nodes );

    // For those nodes that are in the box, get the elements that they
    // construct. These elements are in the box.
    ordinal_type num_elements = std::distance( MT::elementsBegin( mesh ),
					       MT::elementsEnd( mesh ) );
    std::size_t nodes_per_element = MT::nodesPerElement( mesh );
    ordinal_type node_index;
    handle_type node_handle;
    int this_element_in_box;
    typename MT::const_connectivity_iterator mesh_connectivity =
	MT::connectivityBegin( mesh );
    for ( ordinal_type n = 0; n < num_elements; ++n )
    {
	this_element_in_box = 0;
	for ( std::size_t i = 0; i < nodes_per_element; ++i )
	{
	    node_handle = mesh_connectivity[ i*num_elements + n ];
	    node_index = node_indices.find( node_handle )->second;
	    if ( nodes_in_box[ node_index ] )
	    {
		this_element_in_box = 1;
	    }
	}
	elements_in_box.push_back( this_element_in_box );
    }
    assert( (ordinal_type) elements_in_box.size() == num_elements );

    // Get the nodes that belong to the elements in the box, but are not in
    // the box themselves. These will also be used in RCB.
    for ( ordinal_type n = 0; n < num_elements; ++n )
    {
	if ( elements_in_box[n] )
	{
	    for ( std::size_t i = 0; i < nodes_per_element; ++i )
	    {
		node_handle = mesh_connectivity[ i*num_elements + n ];
		node_index = node_indices.find( node_handle )->second;
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
template<class Mesh>
void Rendezvous<Mesh>::sendMeshToRendezvous( 
    const Mesh& mesh, const Teuchos::Array<int>& elements_in_box )
{
    // Setup the communication patterns for moving the mesh to the rendezvous
    // decomposition. This will also move the node and element handles to the
    // rendezvous decomposition.
    Teuchos::Array<handle_type> rendezvous_nodes;
    Teuchos::Array<handle_type> rendezvous_elements;
    setupImportCommunication( mesh, elements_in_box,
			      rendezvous_nodes, rendezvous_elements );

    // Setup export node map.
    ordinal_type num_nodes = std::distance( MT::nodesBegin( mesh ), 
					    MT::nodesEnd( mesh ) );
    Teuchos::ArrayView<const handle_type> export_node_view(
	&*MT::nodesBegin( mesh ), num_nodes );
    RCP_TpetraMap export_node_map = Tpetra::createNonContigMap<handle_type>( 
	export_node_view, d_comm );
    testPostcondition( export_node_map != Teuchos::null,
		       "Error creating node export map." );

    // Setup import node map.
    Teuchos::ArrayView<const handle_type> rendezvous_nodes_view = 
	rendezvous_nodes();
    RCP_TpetraMap import_node_map = Tpetra::createNonContigMap<handle_type>(
	rendezvous_nodes_view, d_comm );
    testPostcondition( import_node_map != Teuchos::null,
		       "Error creating node import map." );

    // Setup export element map.
    ordinal_type num_elements = std::distance( MT::elementsBegin( mesh ), 
					       MT::elementsEnd( mesh ) );
    Teuchos::ArrayView<const handle_type> export_element_view(
	&*MT::elementsBegin( mesh ), num_elements );
    RCP_TpetraMap export_element_map = Tpetra::createNonContigMap<handle_type>(
	export_element_view, d_comm );
    testPostcondition( export_element_map != Teuchos::null,
		       "Error creating element export map." );

    // Setup import element map.
    Teuchos::ArrayView<const handle_type> rendezvous_elements_view =
	rendezvous_elements();
    RCP_TpetraMap import_element_map = Tpetra::createNonContigMap<handle_type>(
	rendezvous_elements_view, d_comm );
    testPostcondition( import_element_map != Teuchos::null,
		       "Error creating element import map." );

    // Setup importers.
    Tpetra::Import<handle_type> node_importer( export_node_map, 
					       import_node_map );
    Tpetra::Import<handle_type> element_importer( export_element_map, 
						  import_element_map );

    // Move the node coordinates to the rendezvous decomposition.
    ordinal_type num_coords = 3*num_nodes;
    Teuchos::ArrayRCP<double> export_coords_view( 
	(double*) &*MT::coordsBegin( mesh ), 0, num_coords, false );
    Teuchos::RCP< Tpetra::MultiVector<double> > export_coords = 
	createMultiVectorFromView( export_node_map, export_coords_view, 
				   num_nodes, 3 );
    Tpetra::MultiVector<double> import_coords( import_node_map, 3 );
    import_coords.doImport( *export_coords, node_importer, Tpetra::INSERT );

    // Move the element connectivity to the rendezvous decomposition.
    int nodes_per_element = MT::nodesPerElement( mesh );
    ordinal_type num_conn = nodes_per_element * num_elements;
    Teuchos::ArrayRCP<handle_type> export_conn_view( 
	(handle_type*) &*MT::connectivityBegin( mesh ),	0, num_conn, false );
    Teuchos::RCP< Tpetra::MultiVector<handle_type,handle_type> > export_conn = 
	createMultiVectorFromView( export_element_map, export_conn_view, 
				   num_elements, nodes_per_element );
    Tpetra::MultiVector<handle_type,handle_type> import_conn( 
	import_element_map, nodes_per_element );
    import_conn.doImport( *export_conn, element_importer, Tpetra::INSERT );

    // Construct the mesh container from the collected data, effectively
    // wrapping it with mesh traits. 
    MeshContainer<handle_type> mesh_container( 
	Teuchos::arcpFromArray( rendezvous_nodes ), 
	import_coords.get1dView(),
	MT::elementType( mesh ),
	MT::elementTopology( mesh ),
	nodes_per_element,
	Teuchos::arcpFromArray( rendezvous_elements ), 
	import_conn.get1dView() );

    // Build the concrete rendezvous mesh from the mesh container.
    d_rendezvous_mesh = createRendezvousMesh( mesh_container );
    testPostcondition( d_rendezvous_mesh != Teuchos::null,
		       "Error creating rendezvous mesh." );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Setup the import communication patterns.
 */
template<class Mesh>
void Rendezvous<Mesh>::setupImportCommunication( 
    const Mesh& mesh,
    const Teuchos::Array<int>& elements_in_box,
    Teuchos::Array<handle_type>& rendezvous_nodes,
    Teuchos::Array<handle_type>& rendezvous_elements )
{
    // Create a node index map for logarithmic time access to connectivity
    // data. 
    typename MT::const_node_iterator export_node_iterator;
    std::map<handle_type,ordinal_type> node_indices;
    ordinal_type m = 0;
    for ( export_node_iterator = MT::nodesBegin( mesh );
	  export_node_iterator != MT::nodesEnd( mesh );
	  ++export_node_iterator )
    {
	node_indices[ *export_node_iterator ] = m;
	++m;
    }

    // Create a element index map for logarithmic time access to connectivity
    // data. 
    typename MT::const_element_iterator export_element_iterator;
    std::map<handle_type,ordinal_type> element_indices;
    m = 0;
    for ( export_element_iterator = MT::elementsBegin( mesh );
	  export_element_iterator != MT::elementsEnd( mesh );
	  ++export_element_iterator )
    {
	element_indices[ *export_element_iterator ] = m;
	++m;
    }

    // Get destination procs for all local elements in the global bounding
    // box. The element will need to be sent to each partition that its
    // connecting nodes exist in. We'll make a unique destination proc set for
    // each element.
    ordinal_type num_nodes = std::distance( MT::nodesBegin( mesh ),
					    MT::nodesEnd( mesh ) );
    ordinal_type num_elements = std::distance( MT::elementsBegin( mesh ),
					       MT::elementsEnd( mesh ) );
    Teuchos::Array< std::set<int> > export_element_procs_set( num_elements );
    std::size_t nodes_per_element = MT::nodesPerElement( mesh );
    ordinal_type node_index;
    handle_type node_handle;
    int destination_proc;
    double node_coords[3];
    typename MT::const_coordinate_iterator mesh_coords = MT::coordsBegin( mesh );
    typename MT::const_connectivity_iterator mesh_connectivity = 
	MT::connectivityBegin( mesh );
    for ( ordinal_type n = 0; n < num_elements; ++n )
    {
	if ( elements_in_box[n] )
	{
	    for ( std::size_t i = 0; i < nodes_per_element; ++i )
	    {
		node_handle = mesh_connectivity[ i*num_elements + n ];
		node_index = node_indices.find( node_handle )->second;
		node_coords[0] = mesh_coords[ node_index ];
		node_coords[1] = mesh_coords[ num_nodes + node_index ];
		node_coords[2] = mesh_coords[ 2*num_nodes + node_index ];
		destination_proc = d_rcb->getDestinationProc( node_coords );
		export_element_procs_set[n].insert( destination_proc );
	    }
	}
    }

    // Unroll the vector of sets into two vectors; one containing the element
    // handle and the other containing the corresponding element destination.
    Teuchos::Array<handle_type> export_elements;
    Teuchos::Array<int> export_element_procs;
    typename MT::const_element_iterator element_iterator;
    Teuchos::Array< std::set<int> >::const_iterator element_vec_iterator;
    std::set<int>::const_iterator element_proc_set_iterator;
    for ( element_vec_iterator = export_element_procs_set.begin(),
	      element_iterator = MT::elementsBegin( mesh );
	  element_vec_iterator != export_element_procs_set.end();
	  ++element_vec_iterator, ++element_iterator )
    {
	for ( element_proc_set_iterator = element_vec_iterator->begin();
	      element_proc_set_iterator != element_vec_iterator->end();
	      ++element_proc_set_iterator )
	{
	    export_elements.push_back( *element_iterator );
	    export_element_procs.push_back( *element_proc_set_iterator );
	}
    }
    export_element_procs_set.clear();

    // Now we know where the elements need to go. Move the elements to the
    // rendezvous decomposition through an inverse communciation operation.
    Tpetra::Distributor element_distributor( d_comm );
    Teuchos::ArrayView<int> export_element_procs_view = export_element_procs();
    ordinal_type num_import_elements = element_distributor.createFromSends(
	export_element_procs_view );
    Teuchos::ArrayView<const handle_type> export_elements_view =
	export_elements();
    Teuchos::Array<handle_type> import_elements( num_import_elements );
    Teuchos::ArrayView<handle_type> import_elements_view = import_elements();
    element_distributor.doPostsAndWaits( 
	export_elements_view, 1, import_elements_view );
    
    // Next move these into the rendezvous element set so that we have a
    // unique list of the elements.
    typename Teuchos::Array<handle_type>::const_iterator import_element_iterator;
    std::set<handle_type> rendezvous_elements_set;
    for ( import_element_iterator = import_elements.begin();
	  import_element_iterator != import_elements.end();
	  ++import_element_iterator )
    {
	rendezvous_elements_set.insert( *import_element_iterator );
    }
    import_elements.clear();

    // Finally put the elements in a Teuchos::Array and get rid of the set.
    rendezvous_elements.resize( rendezvous_elements_set.size() );
    std::copy( rendezvous_elements_set.begin(), rendezvous_elements_set.end(),
	       rendezvous_elements.begin() );
    rendezvous_elements_set.clear();

    // Now get the destination procs for all the nodes. This will be the same
    // destination procs as all of their parent elements. Therefore, nodes may
    // then also have to go to multiple procs because of this and these procs
    // may be different than their original RCB procs.
    ordinal_type element_index;
    typename Teuchos::Array<handle_type>::const_iterator 
	export_elements_iterator;
    Teuchos::Array<int>::const_iterator export_element_procs_iterator;
    Teuchos::Array< std::set<int> > export_node_procs_set( num_nodes );
    for ( export_elements_iterator = export_elements.begin(),
     export_element_procs_iterator = export_element_procs.begin();
	  export_elements_iterator != export_elements.end();
	  ++export_elements_iterator, ++export_element_procs_iterator )
    {
	element_index = 
	    element_indices.find( *export_elements_iterator )->second;

	for ( std::size_t i = 0; i < nodes_per_element; ++i )
	{
	    node_handle = mesh_connectivity[ i*num_elements + element_index ];
	    node_index = node_indices.find( node_handle )->second;

	    export_node_procs_set[ node_index ].insert( 
		*export_element_procs_iterator );
	}
    }
    export_elements.clear();    
    export_element_procs.clear();
    node_indices.clear();
    element_indices.clear();

    // Unroll the vector of sets into two vectors; one containing the node
    // handle and the other containing the corresponding node destination.
    Teuchos::Array<handle_type> export_nodes;
    Teuchos::Array<int> export_node_procs;
    Teuchos::Array< std::set<int> >::const_iterator node_vec_iterator;
    std::set<int>::const_iterator node_proc_set_iterator;
    for ( node_vec_iterator = export_node_procs_set.begin(),
       export_node_iterator = MT::nodesBegin( mesh );
	  node_vec_iterator != export_node_procs_set.end();
	  ++node_vec_iterator, ++export_node_iterator )
    {
	for ( node_proc_set_iterator = node_vec_iterator->begin();
	      node_proc_set_iterator != node_vec_iterator->end();
	      ++node_proc_set_iterator )
	{
	    export_nodes.push_back( *export_node_iterator );
	    export_node_procs.push_back( *node_proc_set_iterator );
	}
    }
    export_node_procs_set.clear();

    // Now we know where the nodes need to go. Move the nodes to the
    // rendezvous decomposition through an inverse communciation operation.
    Tpetra::Distributor node_distributor( d_comm );
    Teuchos::ArrayView<int> export_node_procs_view = export_node_procs();
    ordinal_type num_import_nodes = node_distributor.createFromSends(
	export_node_procs_view );
    Teuchos::ArrayView<const handle_type> export_nodes_view = export_nodes();
    Teuchos::Array<handle_type> import_nodes( num_import_nodes );
    Teuchos::ArrayView<handle_type> import_nodes_view = import_nodes();
    node_distributor.doPostsAndWaits( export_nodes_view, 1, import_nodes_view );
    export_nodes.clear();
    export_node_procs.clear();

    // Next move these into the rendezvous node set so that we have a unique
    // list of the nodes.
    typename Teuchos::Array<handle_type>::const_iterator import_node_iterator;
    std::set<handle_type> rendezvous_nodes_set;
    for ( import_node_iterator = import_nodes.begin();
	  import_node_iterator != import_nodes.end();
	  ++import_node_iterator )
    {
	rendezvous_nodes_set.insert( *import_node_iterator );
    }
    import_nodes.clear();

    // Finally put the nodes in a Teuchos::Array and get rid of the set.
    rendezvous_nodes.resize( rendezvous_nodes_set.size() );
    std::copy( rendezvous_nodes_set.begin(), rendezvous_nodes_set.end(),
	       rendezvous_nodes.begin() );
    rendezvous_nodes_set.clear();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_RENDEZVOUS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_Rendezvous_def.hpp
//---------------------------------------------------------------------------//
