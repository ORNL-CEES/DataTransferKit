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

#include <set>
#include <algorithm>
#include <cassert>

#include "DTK_MeshTools.hpp"
#include <DTK_Exception.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Tuple.hpp>

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
void Rendezvous<Mesh>::build( const RCP_MeshManager& mesh_manager )
{
    // Get the dimension for the mesh.
    d_node_dim = mesh_manager->dim();

    // Extract the mesh nodes and elements that are in the bounding box. These
    // are the pieces of the mesh that will be repartitioned.
    getMeshInBox( mesh_manager );

    // Construct the rendezvous partitioning for the mesh with RCB using the
    // nodes that are in the box.
    d_rcb = Teuchos::rcp( new RCB<Mesh>( mesh_manager ) );
    testPostcondition( d_rcb != Teuchos::null,
		       "Error creating RCB decomposition." );
    d_rcb->partition();

    // Send the mesh in the box to the rendezvous decomposition and build the
    // concrete mesh blocks.
    MeshManager<MeshContainerType> rendezvous_mesh_manager =
	sendMeshToRendezvous( mesh_manager );

    // Build the concrete rendezvous mesh from the mesh container.
    d_rendezvous_mesh = createRendezvousMesh( rendezvous_mesh_manager );
    testPostcondition( d_rendezvous_mesh != Teuchos::null,
		       "Error creating rendezvous mesh." );

    // Create a kD-tree in the rendezvous decomposition.
    d_kdtree = Teuchos::rcp( 
	new KDTree<GlobalOrdinal>( d_rendezvous_mesh , d_node_dim ) );
    testPostcondition( d_kdtree != Teuchos::null,
		       "Error creating rendezvous kD-tree." );
    d_kdtree->build();
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Get the rendezvous destination processes for a blocked list of
 * coordinates that are in the primary decomposition.
 */
template<class Mesh>
Teuchos::Array<int> Rendezvous<Mesh>::procsContainingPoints(
    const Teuchos::ArrayRCP<double>& coords ) const
{
    Teuchos::Array<double> point( d_node_dim );
    GlobalOrdinal num_points = coords.size() / d_node_dim;
    int rendezvous_proc;
    Teuchos::Array<int> destination_procs( num_points );
    for ( GlobalOrdinal n = 0; n < num_points; ++n )
    {
	for ( std::size_t d = 0; d < d_node_dim; ++d )
	{
	    point[d] = coords[ d*num_points + n ];
	}
	rendezvous_proc = d_rcb->getDestinationProc( point );
	destination_procs[n] = rendezvous_proc;
    }

    return destination_procs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the native mesh elements in the rendezvous decomposition
 * containing a blocked list of coordinates also in the rendezvous
 * decomposition. If a point is not found in an element, return an invalid
 * element ordinal and source decomposition proc, -1, for that point.
 */
template<class Mesh>
void Rendezvous<Mesh>::elementsContainingPoints( 
    const Teuchos::ArrayRCP<double>& coords,
    Teuchos::Array<GlobalOrdinal>& elements,
    Teuchos::Array<int>& element_src_procs ) const
{
    Teuchos::Array<double> point( d_node_dim );
    GlobalOrdinal element_ordinal;
    GlobalOrdinal num_points = coords.size() / d_node_dim;
    bool found_point;
    elements.resize( num_points );
    element_src_procs.resize( num_points );
    for ( GlobalOrdinal n = 0; n < num_points; ++n )
    {
	for ( std::size_t d = 0; d < d_node_dim; ++d )
	{
	    point[d] = coords[ d*num_points + n ];
	}

	found_point = d_kdtree->findPoint( point, element_ordinal );

	if ( found_point )
	{
	    elements[n] = element_ordinal;
	    element_src_procs[n] = 
		d_element_src_procs_map.find( element_ordinal )->second;
	}
	else
	{
	    elements[n] = -1;
	    element_src_procs[n] = -1;
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Extract the mesh nodes and elements that are in a bounding box.
 */
template<class Mesh>
void Rendezvous<Mesh>::getMeshInBox( const RCP_MeshManager& mesh_manager )
{
    // Expand the box by a typical mesh element length in all directions plus
    // some tolerance. Doing this catches a few corner cases.
    double tol = 1.0e-4;
    Teuchos::Tuple<double,6> box_bounds = d_global_box.getBounds();
    double box_volume = d_global_box.volume( d_node_dim );
    GlobalOrdinal global_num_elements = mesh_manager->globalNumElements();
    double pow = 1 / d_node_dim;
    double typical_length = std::pow( box_volume / global_num_elements, pow );
    for ( std::size_t d = 0; d < d_node_dim; ++d )
    {
	box_bounds[d] -= typical_length + tol;
	box_bounds[d+3] += typical_length + tol;
    }
    d_global_box = BoundingBox( box_bounds );

    // For every mesh block, get its nodes and elements that are in the
    // expanded box.
    Teuchos::Array<short int> nodes_in_box, elements_in_box;
    BlockIterator block_iterator;
    for ( block_iterator = mesh_manager->blocksBegin();
	  block_iterator != mesh_manager->blocksEnd();
	  ++block_iterator )
    {
	// Setup.
	nodes_in_box.clear();
	elements_in_box.clear();
	int block_id = std::distance( mesh_manager->blocksBegin(), 
				      block_iterator );

	// Create a map indexed by node global ordinal containing the actual
	// node ordinal location in the array. This will give us logarithmic
	// time access to connectivity. I should write a more general hash
	// table to improve this access time as I'm using this strategy for
	// most mesh operations.
	GlobalOrdinal num_nodes = MeshTools<Mesh>::numNodes( *block_iterator );
	std::map<GlobalOrdinal,GlobalOrdinal> node_indices;
	typename MT::const_node_iterator node_iterator;
	GlobalOrdinal array_index = 0;
	for ( node_iterator = MT::nodesBegin( *block_iterator );
	      node_iterator != MT::nodesEnd( *block_iterator );
	      ++node_iterator )
	{
	    node_indices[ *node_iterator ] = array_index;
	    ++array_index;
	}

	// Get all of the nodes that are in the box. 
	Teuchos::Array<double> node_coords( d_node_dim );
	Teuchos::ArrayRCP<const double> mesh_coords =
	    MeshTools<Mesh>::coordsView( *block_iterator );
	for ( GlobalOrdinal n = 0; n < num_nodes; ++n )
	{
	    for ( std::size_t d = 0; d < d_node_dim; ++d )
	    {
		node_coords[d] = mesh_coords[ d*num_nodes + n ];
	    }
	    nodes_in_box.push_back( d_global_box.pointInBox( node_coords ) );
	}
	assert( (GlobalOrdinal) nodes_in_box.size() == num_nodes );

	// For those nodes that are in the box, get the elements that they
	// construct. These elements are in the box.
	GlobalOrdinal num_elements = 
	    MeshTools<Mesh>::numElements( *block_iterator );
	std::size_t nodes_per_element = MT::nodesPerElement( *block_iterator );
	GlobalOrdinal node_index;
	GlobalOrdinal node_ordinal;
	int this_element_in_box;
	Teuchos::ArrayRCP<const GlobalOrdinal> mesh_connectivity = 
	    MeshTools<Mesh>::connectivityView( *block_iterator );
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

	// Get the nodes that belong to the elements in the box, but are not
	// necessarily in the box themselves. These will also be used in RCB.
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

	// Set the active node/element data in the manager for the block.
	mesh_manager->setActiveNodes( nodes_in_box, block_id );
	mesh_manager->setActiveElements( elements_in_box, block_id );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send the mesh to the rendezvous decomposition and rebuild the mesh
 * blocks.
 */
template<class Mesh>
MeshManager<typename Rendezvous<Mesh>::MeshContainerType> 
Rendezvous<Mesh>::sendMeshToRendezvous( 
    const RCP_MeshManager& mesh_manager )
{
    // Setup the mesh blocks.
    Teuchos::ArrayRCP<MeshContainerType> 
	block_containers( mesh_manager->getNumBlocks() );
    BlockIterator block_iterator;
    for ( block_iterator = mesh_manager->blocksBegin();
	  block_iterator != mesh_manager->blocksEnd();
	  ++block_iterator )
    {
	int block_id = std::distance( mesh_manager->blocksBegin(), 
				      block_iterator );

	// Setup the communication patterns for moving the mesh block to the
	// rendezvous decomposition. This will also move the node and element
	// global ordinals to the rendezvous decomposition.
	Teuchos::Array<GlobalOrdinal> rendezvous_nodes;
	Teuchos::Array<GlobalOrdinal> rendezvous_elements;
	setupImportCommunication( *block_iterator, 
				  mesh_manager->getActiveElements( block_id ),
				  rendezvous_nodes, rendezvous_elements );

	// Setup export node map.
	GlobalOrdinal num_nodes = MeshTools<Mesh>::numNodes( *block_iterator );
	Teuchos::ArrayRCP<const GlobalOrdinal> export_node_arcp =
	    MeshTools<Mesh>::nodesView( *block_iterator );
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
	GlobalOrdinal num_elements = 
	    MeshTools<Mesh>::numElements( *block_iterator );
	Teuchos::ArrayRCP<const GlobalOrdinal> export_element_arcp =
	    MeshTools<Mesh>::elementsView( *block_iterator );
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
	Teuchos::ArrayRCP<double> export_coords_view =
	    MeshTools<Mesh>::coordsNonConstView( *block_iterator );
	Teuchos::RCP< Tpetra::MultiVector<double,GlobalOrdinal> > export_coords =
	    Tpetra::createMultiVectorFromView( 
		export_node_map, export_coords_view, num_nodes, d_node_dim );
	Tpetra::MultiVector<double,GlobalOrdinal> 
	    import_coords( import_node_map, d_node_dim );
	import_coords.doImport( *export_coords, node_importer, Tpetra::INSERT );

	// Move the element connectivity to the rendezvous decomposition.
	int nodes_per_element = MT::nodesPerElement( *block_iterator );
	GlobalOrdinal num_conn = nodes_per_element * num_elements;
	Teuchos::ArrayRCP<GlobalOrdinal> export_conn_view =
	    MeshTools<Mesh>::connectivityNonConstView( *block_iterator );
	Teuchos::RCP< Tpetra::MultiVector<GlobalOrdinal,GlobalOrdinal> > 
	    export_conn 
	    = Tpetra::createMultiVectorFromView( 
		export_element_map, export_conn_view, 
		num_elements, nodes_per_element );
	Tpetra::MultiVector<GlobalOrdinal,GlobalOrdinal> import_conn( 
	    import_element_map, nodes_per_element );
	import_conn.doImport( *export_conn, element_importer, Tpetra::INSERT );

	// Construct the mesh block container from the collected data,
	// effectively wrapping it with mesh traits.
	Teuchos::ArrayRCP<GlobalOrdinal> 
	    rendezvous_nodes_array( rendezvous_nodes.size() );
	std::copy( rendezvous_nodes.begin(), rendezvous_nodes.end(),
		   rendezvous_nodes_array.begin() );
	rendezvous_nodes.clear();

	Teuchos::ArrayRCP<GlobalOrdinal> 
	    rendezvous_elements_array( rendezvous_elements.size() );
	std::copy( rendezvous_elements.begin(), rendezvous_elements.end(),
		   rendezvous_elements_array.begin() );
	rendezvous_elements.clear();

	Teuchos::ArrayRCP<const std::size_t> permutation_list = 
	    MeshTools<Mesh>::permutationView( *block_iterator );

	block_containers[ block_id ] = 
	    MeshContainerType( d_node_dim,
			       rendezvous_nodes_array,
			       import_coords.get1dView(),
			       MT::elementTopology( *block_iterator ),
			       nodes_per_element,
			       rendezvous_elements_array,
			       import_conn.get1dView(),
			       permutation_list );
    }

    // Build the rendezvous mesh manager from the rendezvous mesh blocks.
    return MeshManager<MeshContainerType>( block_containers,
					   d_comm,
					   d_node_dim );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Setup the import communication patterns for moving mesh from the
 * primary decomposition to the rendezvous decomposition.
 */
template<class Mesh>
void Rendezvous<Mesh>::setupImportCommunication( 
    const Mesh& mesh,
    const Teuchos::ArrayView<short int>& elements_in_box,
    Teuchos::Array<GlobalOrdinal>& rendezvous_nodes,
    Teuchos::Array<GlobalOrdinal>& rendezvous_elements )
{
    // Create a node index map for logarithmic time access to connectivity
    // data. 
    typename MT::const_node_iterator export_node_iterator;
    std::map<GlobalOrdinal,GlobalOrdinal> node_indices;
    GlobalOrdinal array_index = 0;
    for ( export_node_iterator = MT::nodesBegin( mesh );
	  export_node_iterator != MT::nodesEnd( mesh );
	  ++export_node_iterator )
    {
	node_indices[ *export_node_iterator ] = array_index;
	++array_index;
    }

    // Create a element index map for logarithmic time access to connectivity
    // data. 
    typename MT::const_element_iterator export_element_iterator;
    std::map<GlobalOrdinal,GlobalOrdinal> element_indices;
    array_index = 0;
    for ( export_element_iterator = MT::elementsBegin( mesh );
	  export_element_iterator != MT::elementsEnd( mesh );
	  ++export_element_iterator )
    {
	element_indices[ *export_element_iterator ] = array_index;
	++array_index;
    }

    // Get destination procs for all local elements in the global bounding
    // box. The element will need to be sent to each partition that its
    // connecting nodes exist in. We'll make a unique destination proc set for
    // each element.
    GlobalOrdinal num_nodes = MeshTools<Mesh>::numNodes( mesh );
    GlobalOrdinal num_elements = MeshTools<Mesh>::numElements( mesh );
    Teuchos::Array< std::set<int> > export_element_procs_set( num_elements );
    std::size_t nodes_per_element = MT::nodesPerElement( mesh );
    GlobalOrdinal node_index;
    GlobalOrdinal node_ordinal;
    int destination_proc;
    Teuchos::Array<double> node_coords( d_node_dim );
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
    element_distributor.doPostsAndWaits( export_elements_view, 1, 
					 import_elements_view );
    
    // Extract the rendezvous element source procs from the distributor.
    Teuchos::ArrayView<const int> from_images = element_distributor.getImagesFrom();
    Teuchos::ArrayView<const std::size_t> from_lengths = 
	element_distributor.getLengthsFrom();
    Teuchos::Array<int> element_src_procs;
    for ( int i = 0; i < (int) from_images.size(); ++i )
    {
	for ( std::size_t j = 0; j < from_lengths[i]; ++j )
	{
	    element_src_procs.push_back( from_images[i] );
	}
    }
    testInvariant( element_src_procs.size() == num_import_elements,
		   "number of element src procs != number of import elements" );
        
    // Next, move these into the rendezvous element set so that we have a
    // unique list of the elements and build the rendezvous mesh element to
    // source proc map.
    std::set<GlobalOrdinal> rendezvous_elements_set;
    for ( int n = 0; n < num_import_elements; ++n )
    {
	if ( rendezvous_elements_set.insert( import_elements[n] ).second )
	{
	    d_element_src_procs_map[ import_elements[n] ] = 
		element_src_procs[n];
	}
    }
    import_elements.clear();
    element_src_procs.clear();

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
