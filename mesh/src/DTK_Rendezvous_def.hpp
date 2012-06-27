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
#include "DTK_MeshTools.hpp"
#include <DTK_Exception.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_ArrayView.hpp>

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
    , d_node_dim( 0 )
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
    // Get the node dimension for the mesh.
    d_node_dim = MT::nodeDim( mesh );

    // Extract the mesh nodes and elements that are in the bounding box.
    Teuchos::Array<int> nodes_in_box;
    Teuchos::Array<int> elements_in_box;
    getMeshInBox( mesh, d_global_box, nodes_in_box, elements_in_box );

    // Construct the rendezvous decomposition of the mesh with RCB using the
    // nodes that are in the box.
    d_rcb = Teuchos::rcp( 
	new RCB<Mesh>( mesh, Teuchos::arcpFromArray( nodes_in_box ), 
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
    d_kdtree = Teuchos::rcp( new KDTree<GlobalOrdinal>( d_rendezvous_mesh ) );
    testPostcondition( d_kdtree != Teuchos::null,
		       "Error creating rendezvous kD-tree." );
    d_kdtree->build();
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Get the rendezvous processes for a blocked list of coordinates.
 */
template<class Mesh>
Teuchos::Array<int> Rendezvous<Mesh>::getRendezvousProcs(
    const Teuchos::ArrayRCP<double>& coords ) const
{
    double point[3];
    GlobalOrdinal num_points = coords.size() / d_node_dim;
    int rendezvous_proc;
    Teuchos::Array<int> destination_procs( num_points );
    for ( GlobalOrdinal n = 0; n < num_points; ++n )
    {
	for ( std::size_t d = 0; d < d_node_dim; ++d )
	{
	    point[d] = coords[ d*num_points + n ];
	}
	for ( std::size_t d = d_node_dim; d < 3; ++d )
	{
	    point[d] = 0.0;
	}
	rendezvous_proc = d_rcb->getDestinationProc( point );
	destination_procs[n] = rendezvous_proc;
    }

    testPostcondition( static_cast<GlobalOrdinal>( destination_procs.size() )
		       == num_points,
		       "Error getting destination processes." );

    return destination_procs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the native mesh elements containing a blocked list of
 * coordinates.
 */
template<class Mesh>
Teuchos::Array<typename Rendezvous<Mesh>::GlobalOrdinal>
Rendezvous<Mesh>::getElements( const Teuchos::ArrayRCP<double>& coords ) const
{
    Teuchos::Array<double> point( d_node_dim );
    GlobalOrdinal element_ordinal;
    GlobalOrdinal num_points = coords.size() / d_node_dim;
    Teuchos::Array<GlobalOrdinal> element_ordinals( num_points );
    for ( GlobalOrdinal n = 0; n < num_points; ++n )
    {
	for ( std::size_t d = 0; d < d_node_dim; ++d )
	{
	    point[d] = coords[ d*num_points + n ];
	}
	element_ordinal = d_kdtree->findPoint( point );
	element_ordinals[n] = element_ordinal;
    }

    testPostcondition( static_cast<GlobalOrdinal>( element_ordinals.size() )
		       == num_points,
		       "Error getting mesh elements." );

    return element_ordinals;
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
    // Create a map indexed by node global ordinal containing the actual node
    // ordinal location. This will give us logarithmic time access to
    // connectivity. I should write a more general hash table to improve this
    // access time as I'm using this strategy for most mesh operations.
    GlobalOrdinal num_nodes = std::distance( MT::nodesBegin( mesh ),
					     MT::nodesEnd( mesh ) );
    std::map<GlobalOrdinal,GlobalOrdinal> node_indices;
    typename MT::const_node_iterator node_iterator;
    GlobalOrdinal m = 0;
    for ( node_iterator = MT::nodesBegin( mesh );
	  node_iterator != MT::nodesEnd( mesh );
	  ++node_iterator )
    {
	node_indices[ *node_iterator ] = m;
	++m;
    }

    // Get all of the nodes that are in the box. 
    double node_coords[3];
    Teuchos::ArrayRCP<const double> mesh_coords =
	MeshTools<Mesh>::coordsView( mesh );
    for ( GlobalOrdinal n = 0; n < num_nodes; ++n )
    {
	for ( std::size_t d = 0; d < d_node_dim; ++d )
	{
	    node_coords[d] = mesh_coords[ d*num_nodes + n ];
	}
	for ( std::size_t d = d_node_dim; d < 3; ++d )
	{
	    node_coords[d] = 0.0;
	}
	nodes_in_box.push_back( box.pointInBox( node_coords ) );
    }
    assert( (GlobalOrdinal) nodes_in_box.size() == num_nodes );

    // For those nodes that are in the box, get the elements that they
    // construct. These elements are in the box.
    GlobalOrdinal num_elements = std::distance( MT::elementsBegin( mesh ),
						MT::elementsEnd( mesh ) );
    std::size_t nodes_per_element = MT::nodesPerElement( mesh );
    GlobalOrdinal node_index;
    GlobalOrdinal node_ordinal;
    int this_element_in_box;
    Teuchos::ArrayRCP<const GlobalOrdinal> mesh_connectivity = 
	MeshTools<Mesh>::connectivityView( mesh );
    for ( GlobalOrdinal n = 0; n < num_elements; ++n )
    {
	this_element_in_box = 0;
	for ( std::size_t i = 0; i < nodes_per_element; ++i )
	{
	    node_ordinal = mesh_connectivity[ i*num_elements + n ];
	    node_index = node_indices.find( node_ordinal )->second;
	    if ( nodes_in_box[ node_index ] )
	    {
		this_element_in_box = 1;
	    }
	}
	elements_in_box.push_back( this_element_in_box );
    }
    assert( (GlobalOrdinal) elements_in_box.size() == num_elements );

    // Get the nodes that belong to the elements in the box, but are not in
    // the box themselves. These will also be used in RCB.
    for ( GlobalOrdinal n = 0; n < num_elements; ++n )
    {
	if ( elements_in_box[n] )
	{
	    for ( std::size_t i = 0; i < nodes_per_element; ++i )
	    {
		node_ordinal = mesh_connectivity[ i*num_elements + n ];
		node_index = node_indices.find( node_ordinal )->second;
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
    // decomposition. This will also move the node and element global ordinals
    // to the rendezvous decomposition.
    Teuchos::Array<GlobalOrdinal> rendezvous_nodes;
    Teuchos::Array<GlobalOrdinal> rendezvous_elements;
    setupImportCommunication( mesh, elements_in_box,
			      rendezvous_nodes, rendezvous_elements );

    // Setup export node map.
    GlobalOrdinal num_nodes = std::distance( MT::nodesBegin( mesh ), 
					     MT::nodesEnd( mesh ) );

    Teuchos::ArrayRCP<const GlobalOrdinal> export_node_arcp =
	MeshTools<Mesh>::nodesView( mesh );
    Teuchos::ArrayView<const GlobalOrdinal> export_node_view =
	export_node_arcp();
    RCP_TpetraMap export_node_map = Tpetra::createNonContigMap<GlobalOrdinal>( 
	export_node_view, d_comm );
    testPostcondition( export_node_map != Teuchos::null,
		       "Error creating node export map." );

    // Setup import node map.
    Teuchos::ArrayView<const GlobalOrdinal> rendezvous_nodes_view = 
	rendezvous_nodes();
    RCP_TpetraMap import_node_map = Tpetra::createNonContigMap<GlobalOrdinal>(
	rendezvous_nodes_view, d_comm );
    testPostcondition( import_node_map != Teuchos::null,
		       "Error creating node import map." );

    // Setup export element map.
    GlobalOrdinal num_elements = std::distance( MT::elementsBegin( mesh ), 
						MT::elementsEnd( mesh ) );
    Teuchos::ArrayRCP<const GlobalOrdinal> export_element_arcp =
	MeshTools<Mesh>::elementsView( mesh );
    Teuchos::ArrayView<const GlobalOrdinal> export_element_view =
	export_element_arcp();
    RCP_TpetraMap export_element_map = 
	Tpetra::createNonContigMap<GlobalOrdinal>(
	    export_element_view, d_comm );
    testPostcondition( export_element_map != Teuchos::null,
		       "Error creating element export map." );

    // Setup import element map.
    Teuchos::ArrayView<const GlobalOrdinal> rendezvous_elements_view =
	rendezvous_elements();
    RCP_TpetraMap import_element_map = 
	Tpetra::createNonContigMap<GlobalOrdinal>(
	    rendezvous_elements_view, d_comm );
    testPostcondition( import_element_map != Teuchos::null,
		       "Error creating element import map." );

    // Setup importers.
    Tpetra::Import<GlobalOrdinal> node_importer( export_node_map, 
						 import_node_map );
    Tpetra::Import<GlobalOrdinal> element_importer( export_element_map, 
						    import_element_map );

    // Move the node coordinates to the rendezvous decomposition.
    GlobalOrdinal num_coords = d_node_dim*num_nodes;
    Teuchos::ArrayRCP<double> export_coords_view;
    if ( num_coords == 0 )
    {
	export_coords_view = Teuchos::ArrayRCP<double>( 0, 0.0 );
    }
    else
    {
	export_coords_view = MeshTools<Mesh>::coordsNonConstView( mesh );
    }
    Teuchos::RCP< Tpetra::MultiVector<double,GlobalOrdinal> > export_coords = 
	createMultiVectorFromView( export_node_map, export_coords_view, 
				   num_nodes, d_node_dim );
    Tpetra::MultiVector<double,GlobalOrdinal> 
	import_coords( import_node_map, d_node_dim );
    import_coords.doImport( *export_coords, node_importer, Tpetra::INSERT );

    // Move the element connectivity to the rendezvous decomposition.
    int nodes_per_element = MT::nodesPerElement( mesh );
    GlobalOrdinal num_conn = nodes_per_element * num_elements;
    Teuchos::ArrayRCP<GlobalOrdinal> export_conn_view;
    if ( num_conn == 0 )
    {
	export_conn_view = Teuchos::ArrayRCP<GlobalOrdinal>( 0, 0.0 );
    }
    else
    {
	export_conn_view = MeshTools<Mesh>::connectivityNonConstView( mesh );
    }
    Teuchos::RCP< Tpetra::MultiVector<GlobalOrdinal,GlobalOrdinal> > export_conn 
	= createMultiVectorFromView( export_element_map, export_conn_view, 
				     num_elements, nodes_per_element );
    Tpetra::MultiVector<GlobalOrdinal,GlobalOrdinal> import_conn( 
	import_element_map, nodes_per_element );
    import_conn.doImport( *export_conn, element_importer, Tpetra::INSERT );

    // Construct the mesh container from the collected data, effectively
    // wrapping it with mesh traits.
    Teuchos::ArrayRCP<const std::size_t> permutation_list = 
	MeshTools<Mesh>::permutationView( mesh );
    MeshContainer<GlobalOrdinal> mesh_container( 
	d_node_dim,
	Teuchos::arcpFromArray( rendezvous_nodes ), 
	import_coords.get1dView(),
	MT::elementTopology( mesh ),
	nodes_per_element,
	Teuchos::arcpFromArray( rendezvous_elements ), 
	import_conn.get1dView(),
	permutation_list );

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
    Teuchos::Array<GlobalOrdinal>& rendezvous_nodes,
    Teuchos::Array<GlobalOrdinal>& rendezvous_elements )
{
    // Create a node index map for logarithmic time access to connectivity
    // data. 
    typename MT::const_node_iterator export_node_iterator;
    std::map<GlobalOrdinal,GlobalOrdinal> node_indices;
    GlobalOrdinal m = 0;
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
    std::map<GlobalOrdinal,GlobalOrdinal> element_indices;
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
    GlobalOrdinal num_nodes = std::distance( MT::nodesBegin( mesh ),
					     MT::nodesEnd( mesh ) );
    GlobalOrdinal num_elements = std::distance( MT::elementsBegin( mesh ),
						MT::elementsEnd( mesh ) );
    Teuchos::Array< std::set<int> > export_element_procs_set( num_elements );
    std::size_t nodes_per_element = MT::nodesPerElement( mesh );
    GlobalOrdinal node_index;
    GlobalOrdinal node_ordinal;
    int destination_proc;
    double node_coords[3];
    Teuchos::ArrayRCP<const double> mesh_coords =
	MeshTools<Mesh>::coordsView( mesh );
    Teuchos::ArrayRCP<const GlobalOrdinal> mesh_connectivity =
	MeshTools<Mesh>::connectivityView( mesh );
    for ( GlobalOrdinal n = 0; n < num_elements; ++n )
    {
	if ( elements_in_box[n] )
	{
	    for ( std::size_t i = 0; i < nodes_per_element; ++i )
	    {
		node_ordinal = mesh_connectivity[ i*num_elements + n ];
		node_index = node_indices.find( node_ordinal )->second;
		for ( std::size_t d = 0; d < d_node_dim; ++d )
		{
		    node_coords[d] = mesh_coords[ d*num_nodes + node_index ];
		}
		for ( std::size_t d = d_node_dim; d < 3; ++d )
		{
		    node_coords[d] = 0.0;
		}
		destination_proc = d_rcb->getDestinationProc( node_coords );
		export_element_procs_set[n].insert( destination_proc );
	    }
	}
    }

    // Unroll the vector of sets into two vectors; one containing the element
    // ordinal and the other containing the corresponding element destination.
    Teuchos::Array<GlobalOrdinal> export_elements;
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
    GlobalOrdinal num_import_elements = element_distributor.createFromSends(
	export_element_procs_view );
    Teuchos::ArrayView<const GlobalOrdinal> export_elements_view =
	export_elements();
    Teuchos::Array<GlobalOrdinal> import_elements( num_import_elements );
    Teuchos::ArrayView<GlobalOrdinal> import_elements_view = import_elements();
    element_distributor.doPostsAndWaits( 
	export_elements_view, 1, import_elements_view );
    
    // Next move these into the rendezvous element set so that we have a
    // unique list of the elements.
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
	import_element_iterator;
    std::set<GlobalOrdinal> rendezvous_elements_set;
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
    GlobalOrdinal element_index;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
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
	    node_ordinal = mesh_connectivity[ i*num_elements + element_index ];
	    node_index = node_indices.find( node_ordinal )->second;

	    export_node_procs_set[ node_index ].insert( 
		*export_element_procs_iterator );
	}
    }
    export_elements.clear();    
    export_element_procs.clear();
    node_indices.clear();
    element_indices.clear();

    // Unroll the vector of sets into two vectors; one containing the node
    // ordinal and the other containing the corresponding node destination.
    Teuchos::Array<GlobalOrdinal> export_nodes;
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
    GlobalOrdinal num_import_nodes = node_distributor.createFromSends(
	export_node_procs_view );
    Teuchos::ArrayView<const GlobalOrdinal> export_nodes_view = export_nodes();
    Teuchos::Array<GlobalOrdinal> import_nodes( num_import_nodes );
    Teuchos::ArrayView<GlobalOrdinal> import_nodes_view = import_nodes();
    node_distributor.doPostsAndWaits( export_nodes_view, 1, 
				      import_nodes_view );
    export_nodes.clear();
    export_node_procs.clear();

    // Next move these into the rendezvous node set so that we have a unique
    // list of the nodes.
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
	import_node_iterator;
    std::set<GlobalOrdinal> rendezvous_nodes_set;
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
