//---------------------------------------------------------------------------//
/*!
 * \file DTK_Rendezvous_def.hpp
 * \author Stuart R. Slattery
 * \brief Rendezvous definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_RENDEZVOUS_DEF_HPP
#define DTK_RENDEZVOUS_DEF_HPP

#include <set>
#include <map>
#include <algorithm>
#include <cassert>

#include "DTK_Rendezvous.hpp"
#include <DTK_Exception.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Hashtable.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include <Tpetra_Map.hpp>
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
    sendMeshToRendezvous( nodes_in_box, elements_in_box );

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
    const Mesh& mesh,
    const std::vector<char>& nodes_in_box, 
    const std::vector<char>& elements_in_box )
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

    // Move nodes to rendezvous decomposition. This includes coordinates that
    // aren't specified initially by RCB, but are needed to complete the
    // elements that belong to nodes on the RCB boundaries.
    int num_nodes = std::distance( MT::nodesBegin( mesh ), 
				   MT::nodesEnd( mesh ) );
    Teuchos::ArrayView<zoltan_id_type> export_node_ids = 
	d_rcb->getExportGlobalIds();
    Teuchos::ArrayView<zoltan_id_type> export_node_indices =
	d_rcb->getExportLocalIds();
    std::vector<zoltan_id_type> export_node_map_ids( 
	num_nodes, Teuchos::OrdinalTraits<zoltan_id_type>::max() );
    Teuchos::ArrayView<zoltan_id_type>::const_iterator 
	export_node_id_iterator;
    Teuchos::ArrayView<zoltan_id_type>::const_iterator 
	export_node_index_iterator;
    for ( export_node_id_iterator = export_node_ids.begin(),
       export_node_index_iterator = export_node_indices.begin();
	  export_node_id_iterator != export_node_ids.end();
	  ++export_node_id_iterator, ++export_node_index_iterator )
    {
	export_node_map_ids[ *export_node_index_iterator ] =
	    *export_node_ids_iterator;
    }
    Teuchos::RCP< Tpetra::Map<unsiged int> > export_node_map = 
	Tpetra::createNonContigMap<zoltan_id_type>( 
	    Teuchos::ArrayView( export_node_map_ids ), d_global_comm );

    Teuchos::ArrayView<zoltan_id_type> import_node_ids = 
	d_rcb->getImportGlobalIds();
    Teuchos::RCP< Tpetra::Map<zoltan_id_type> > import_node_map = 
	Tpetra::createNonContigMap<zoltan_id_type>( 
	    import_node_ids, d_global_comm );

    Tpetra::Export node_exporter( export_node_map, import_node_map );

    MT::const_handle_iterator export_nodes = MT::nodesBegin( mesh );
    Teuchos::ArrayRCP<handle_type> export_node_view( 
	&*export_nodes, 0, d_rcb->numExport(), false );
    Teuchos::RCP< Tpetra::Vector<handle_type> > export_node_vector =
	Tpetra::createVectorFromView( export_node_map, export_node_view );

    Teuchos::RCP< Tpetra::Vector<handle_type> > import_node_vector =
	Tpetra::createVector( import_node_map );

    import_node_vector->doExport( 
	*export_node_vector, node_exporter, Tpetra::Insert );

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
    export_node_map_ids.clear();

    // Move the elements to the rendezvous decomposition.
    

    // Move the connectivity to the rendezvous decomposition.

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

    // Build the rendezous mesh from the mesh container.
    d_rendezvous_mesh = createRendezvousMesh( mesh_container );

    testPostcondition( d_rendezvous_mesh != Teuchos::null,
		       "Error creating rendezvous mesh." );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_RENDEZVOUS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_Rendezvous_def.hpp
//---------------------------------------------------------------------------//
