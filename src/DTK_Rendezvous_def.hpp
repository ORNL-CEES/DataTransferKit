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

#include "DTK_MeshTools.hpp"
#include "DTK_Assertion.hpp"

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
 *
 * \param comm The communicator over which to build the rendezvous
 * decomposition. 
 *
 * \param global_box The global bounding box inside of which the rendezvous
 * decomposition will be generated.
 */
template<class Mesh>
Rendezvous<Mesh>::Rendezvous( const RCP_Comm& comm,
			      const BoundingBox& global_box )
    : d_comm( comm )
    , d_global_box( global_box )
    , d_vertex_dim( 0 )
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
 *
 * \param mesh_manager The mesh to repartition to the rendezvous
 * decomposition. 
 */
template<class Mesh> 
void Rendezvous<Mesh>::build( const RCP_MeshManager& mesh_manager )
{
    // Get the dimension for the mesh.
    d_vertex_dim = mesh_manager->dim();

    // Extract the mesh vertices and elements that are in the bounding
    // box. These are the pieces of the mesh that will be repartitioned.
    getMeshInBox( mesh_manager );

    // Construct the rendezvous partitioning for the mesh with RCB using the
    // vertices that are in the box.
    d_rcb = Teuchos::rcp( new RCB<Mesh>( mesh_manager ) );
    testPostcondition( d_rcb != Teuchos::null );
    d_rcb->partition();

    // Send the mesh in the box to the rendezvous decomposition and build the
    // concrete mesh blocks.
    MeshManager<MeshContainerType> rendezvous_mesh_manager =
	sendMeshToRendezvous( mesh_manager );

    // Build the concrete rendezvous mesh from the mesh container.
    d_rendezvous_mesh = createRendezvousMesh( rendezvous_mesh_manager );
    testPostcondition( d_rendezvous_mesh != Teuchos::null );

    // Create a kD-tree in the rendezvous decomposition.
    d_kdtree = Teuchos::rcp( 
	new KDTree<GlobalOrdinal>( d_rendezvous_mesh , d_vertex_dim ) );
    testPostcondition( d_kdtree != Teuchos::null );
    d_kdtree->build();
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Get the rendezvous destination processes for a blocked list of
 * coordinates that are in the primary decomposition.
 *
 * \param coords A blocked list of coordinates for which rendezvous
 * decomposition destination procs are desired.
 *
 * \return An array of the rendezvous decomposition destination procs. A proc
 * will be returned for each point in the same order as the points were
 * provided. 
 */
template<class Mesh>
Teuchos::Array<int> Rendezvous<Mesh>::procsContainingPoints(
    const Teuchos::ArrayRCP<double>& coords ) const
{
    Teuchos::Array<double> point( d_vertex_dim );
    GlobalOrdinal num_points = coords.size() / d_vertex_dim;
    int rendezvous_proc;
    Teuchos::Array<int> destination_procs( num_points );
    for ( GlobalOrdinal n = 0; n < num_points; ++n )
    {
	for ( int d = 0; d < d_vertex_dim; ++d )
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
 * decomposition. 
 * 
 * \param coords A blocked list of coordinates to search the mesh with. 

 * \param elements An array of the elements the points were found in. An
 * element will be returned for each point in the order they were provided
 * in. If a point is not found in an element, return an invalid element
 * ordinal, -1, for that point.
 *
 * \param element_src_procs The source procs that own the elements. Once proc
 * is provided for each element in the order that the elements were
 * provided. If a point is not found in an element, return an invalid element
 * source proc, -1, for that point.
 */
template<class Mesh>
void Rendezvous<Mesh>::elementsContainingPoints( 
    const Teuchos::ArrayRCP<double>& coords,
    Teuchos::Array<GlobalOrdinal>& elements,
    Teuchos::Array<int>& element_src_procs ) const
{
    Teuchos::Array<double> point( d_vertex_dim );
    GlobalOrdinal element_ordinal;
    GlobalOrdinal num_points = coords.size() / d_vertex_dim;
    bool found_point;
    elements.resize( num_points );
    element_src_procs.resize( num_points );
    for ( GlobalOrdinal n = 0; n < num_points; ++n )
    {
	for ( int d = 0; d < d_vertex_dim; ++d )
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
 * \brief Extract the mesh vertices and elements that are in a bounding box.
 *
 * \param mesh_manager The mesh to search the box with.
 */
template<class Mesh>
void Rendezvous<Mesh>::getMeshInBox( const RCP_MeshManager& mesh_manager )
{
    // Expand the box by a typical mesh element length in all directions plus
    // some tolerance. Doing this catches a few corner cases.
    double tol = 1.0e-4;
    Teuchos::Tuple<double,6> box_bounds = d_global_box.getBounds();
    double box_volume = d_global_box.volume( d_vertex_dim );
    GlobalOrdinal global_num_elements = mesh_manager->globalNumElements();
    double pow = 1 / d_vertex_dim;
    double typical_length = std::pow( box_volume / global_num_elements, pow );
    for ( int d = 0; d < d_vertex_dim; ++d )
    {
	box_bounds[d] -= typical_length + tol;
	box_bounds[d+3] += typical_length + tol;
    }
    d_global_box = BoundingBox( box_bounds );

    // For every mesh block, get its vertices and elements that are in the
    // expanded box.
    Teuchos::Array<short int> vertices_in_box, elements_in_box;
    BlockIterator block_iterator;
    for ( block_iterator = mesh_manager->blocksBegin();
	  block_iterator != mesh_manager->blocksEnd();
	  ++block_iterator )
    {
	// Setup.
	vertices_in_box.clear();
	elements_in_box.clear();
	int block_id = std::distance( mesh_manager->blocksBegin(), 
				      block_iterator );

	// Create a map indexed by vertex global ordinal containing the actual
	// vertex ordinal location in the array. This will give us logarithmic
	// time access to connectivity. I should write a more general hash
	// table to improve this access time as I'm using this strategy for
	// most mesh operations.
	GlobalOrdinal num_vertices = 
	    MeshTools<Mesh>::numVertices( *block_iterator );
	std::map<GlobalOrdinal,GlobalOrdinal> vertex_indices;
	typename MT::const_vertex_iterator vertex_iterator;
	GlobalOrdinal array_index = 0;
	for ( vertex_iterator = MT::verticesBegin( *block_iterator );
	      vertex_iterator != MT::verticesEnd( *block_iterator );
	      ++vertex_iterator )
	{
	    vertex_indices[ *vertex_iterator ] = array_index;
	    ++array_index;
	}

	// Get all of the vertices that are in the box. 
	Teuchos::Array<double> vertex_coords( d_vertex_dim );
	Teuchos::ArrayRCP<const double> mesh_coords =
	    MeshTools<Mesh>::coordsView( *block_iterator );
	for ( GlobalOrdinal n = 0; n < num_vertices; ++n )
	{
	    for ( int d = 0; d < d_vertex_dim; ++d )
	    {
		vertex_coords[d] = mesh_coords[ d*num_vertices + n ];
	    }
	    vertices_in_box.push_back( 
		d_global_box.pointInBox( vertex_coords ) );
	}
	testInvariant( (GlobalOrdinal) vertices_in_box.size() == num_vertices );

	// For those vertices that are in the box, get the elements that they
	// construct. These elements are in the box.
	GlobalOrdinal num_elements = 
	    MeshTools<Mesh>::numElements( *block_iterator );
	int vertices_per_element = MT::verticesPerElement( *block_iterator );
	GlobalOrdinal vertex_index;
	GlobalOrdinal vertex_ordinal;
	int this_element_in_box;
	Teuchos::ArrayRCP<const GlobalOrdinal> mesh_connectivity = 
	    MeshTools<Mesh>::connectivityView( *block_iterator );
	for ( GlobalOrdinal n = 0; n < num_elements; ++n )
	{
	    this_element_in_box = 0;
	    for ( int i = 0; i < vertices_per_element; ++i )
	    {
		vertex_ordinal = mesh_connectivity[ i*num_elements + n ];
		vertex_index = vertex_indices.find( vertex_ordinal )->second;
		if ( vertices_in_box[ vertex_index ] )
		{
		    this_element_in_box = 1;
		}
	    }
	    elements_in_box.push_back( this_element_in_box );
	}
	testInvariant( (GlobalOrdinal) elements_in_box.size() == num_elements );

	// Get the vertices that belong to the elements in the box, but are not
	// necessarily in the box themselves. These will also be used in RCB.
	for ( GlobalOrdinal n = 0; n < num_elements; ++n )
	{
	    if ( elements_in_box[n] )
	    {
		for ( int i = 0; i < vertices_per_element; ++i )
		{
		    vertex_ordinal = mesh_connectivity[ i*num_elements + n ];
		    vertex_index = 
			vertex_indices.find( vertex_ordinal )->second;
		    vertices_in_box[ vertex_index ] = 1;
		}
	    }
	}

	// Set the active vertex/element data in the manager for the block.
	mesh_manager->setActiveVertices( vertices_in_box, block_id );
	mesh_manager->setActiveElements( elements_in_box, block_id );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send the mesh to the rendezvous decomposition and rebuild the mesh
 * blocks.
 *
 * \param mesh_manager The mesh to send to the rendezvous decomposition.
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
	// rendezvous decomposition. This will also move the vertex and element
	// global ordinals to the rendezvous decomposition.
	Teuchos::Array<GlobalOrdinal> rendezvous_vertices;
	Teuchos::Array<GlobalOrdinal> rendezvous_elements;
	setupImportCommunication( *block_iterator, 
				  mesh_manager->getActiveElements( block_id ),
				  rendezvous_vertices, rendezvous_elements );

	// Setup export vertex map.
	GlobalOrdinal num_vertices = 
	    MeshTools<Mesh>::numVertices( *block_iterator );
	Teuchos::ArrayRCP<const GlobalOrdinal> export_vertex_arcp =
	    MeshTools<Mesh>::verticesView( *block_iterator );
	Teuchos::ArrayView<const GlobalOrdinal> export_vertex_view =
	    export_vertex_arcp();
	RCP_TpetraMap export_vertex_map = 
	    Tpetra::createNonContigMap<GlobalOrdinal>( 
		export_vertex_view, d_comm );
	testInvariant( export_vertex_map != Teuchos::null );

	// Setup import vertex map.
	Teuchos::ArrayView<const GlobalOrdinal> rendezvous_vertices_view = 
	    rendezvous_vertices();
	RCP_TpetraMap import_vertex_map = 
	    Tpetra::createNonContigMap<GlobalOrdinal>(
		rendezvous_vertices_view, d_comm );
	testInvariant( import_vertex_map != Teuchos::null );

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
	testInvariant( export_element_map != Teuchos::null );

	// Setup import element map.
	Teuchos::ArrayView<const GlobalOrdinal> rendezvous_elements_view =
	    rendezvous_elements();
	RCP_TpetraMap import_element_map = 
	    Tpetra::createNonContigMap<GlobalOrdinal>(
		rendezvous_elements_view, d_comm );
	testInvariant( import_element_map != Teuchos::null );

	// Setup importers.
	Tpetra::Import<GlobalOrdinal> vertex_importer( export_vertex_map, 
						       import_vertex_map );
	Tpetra::Import<GlobalOrdinal> element_importer( export_element_map, 
							import_element_map );

	// Move the vertex coordinates to the rendezvous decomposition.
	Teuchos::ArrayRCP<double> export_coords_view =
	    MeshTools<Mesh>::coordsNonConstView( *block_iterator );
	Teuchos::RCP< Tpetra::MultiVector<double,GlobalOrdinal> > 
	    export_coords = Tpetra::createMultiVectorFromView( 
		export_vertex_map, export_coords_view, 
		num_vertices, d_vertex_dim );
	Tpetra::MultiVector<double,GlobalOrdinal> 
	    import_coords( import_vertex_map, d_vertex_dim );
	import_coords.doImport( 
	    *export_coords, vertex_importer, Tpetra::INSERT );

	// Move the element connectivity to the rendezvous decomposition.
	int vertices_per_element = MT::verticesPerElement( *block_iterator );
	Teuchos::ArrayRCP<GlobalOrdinal> export_conn_view =
	    MeshTools<Mesh>::connectivityNonConstView( *block_iterator );
	Teuchos::RCP< Tpetra::MultiVector<GlobalOrdinal,GlobalOrdinal> > 
	    export_conn 
	    = Tpetra::createMultiVectorFromView( 
		export_element_map, export_conn_view, 
		num_elements, vertices_per_element );
	Tpetra::MultiVector<GlobalOrdinal,GlobalOrdinal> import_conn( 
	    import_element_map, vertices_per_element );
	import_conn.doImport( *export_conn, element_importer, Tpetra::INSERT );

	// Construct the mesh block container from the collected data,
	// effectively wrapping it with mesh traits.
	Teuchos::ArrayRCP<GlobalOrdinal> 
	    rendezvous_vertices_array( rendezvous_vertices.size() );
	std::copy( rendezvous_vertices.begin(), rendezvous_vertices.end(),
		   rendezvous_vertices_array.begin() );
	rendezvous_vertices.clear();

	Teuchos::ArrayRCP<GlobalOrdinal> 
	    rendezvous_elements_array( rendezvous_elements.size() );
	std::copy( rendezvous_elements.begin(), rendezvous_elements.end(),
		   rendezvous_elements_array.begin() );
	rendezvous_elements.clear();

	Teuchos::ArrayRCP<const int> permutation_list = 
	    MeshTools<Mesh>::permutationView( *block_iterator );

	block_containers[ block_id ] = 
	    MeshContainerType( d_vertex_dim,
			       rendezvous_vertices_array,
			       import_coords.get1dView(),
			       MT::elementTopology( *block_iterator ),
			       vertices_per_element,
			       rendezvous_elements_array,
			       import_conn.get1dView(),
			       permutation_list );
    }

    // Build the rendezvous mesh manager from the rendezvous mesh blocks.
    return MeshManager<MeshContainerType>( block_containers,
					   d_comm,
					   d_vertex_dim );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Setup the import communication patterns for moving mesh from the
 * primary decomposition to the rendezvous decomposition.
 *
 * \param mesh The mesh block to setup communication for.
 *
 * \param elements_in_box An array of elements in the rendezvous box. Every
 * element gets an entry, 0 if outside the box, 1 if inside.
 *
 * \param rendezvous_vertices An array of the vertex global ordinals in the
 * rendezvous decomposition.
 *
 * \param rendezvous_elements An array of the element global ordinals in the
 * rendezvous decomposition.
 */
template<class Mesh>
void Rendezvous<Mesh>::setupImportCommunication( 
    const Mesh& mesh,
    const Teuchos::ArrayView<short int>& elements_in_box,
    Teuchos::Array<GlobalOrdinal>& rendezvous_vertices,
    Teuchos::Array<GlobalOrdinal>& rendezvous_elements )
{
    // Create a vertex index map for logarithmic time access to connectivity
    // data. 
    typename MT::const_vertex_iterator export_vertex_iterator;
    std::map<GlobalOrdinal,GlobalOrdinal> vertex_indices;
    GlobalOrdinal array_index = 0;
    for ( export_vertex_iterator = MT::verticesBegin( mesh );
	  export_vertex_iterator != MT::verticesEnd( mesh );
	  ++export_vertex_iterator )
    {
	vertex_indices[ *export_vertex_iterator ] = array_index;
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
    // connecting vertices exist in. We'll make a unique destination proc set
    // for each element.
    GlobalOrdinal num_vertices = MeshTools<Mesh>::numVertices( mesh );
    GlobalOrdinal num_elements = MeshTools<Mesh>::numElements( mesh );
    Teuchos::Array< std::set<int> > export_element_procs_set( num_elements );
    int vertices_per_element = MT::verticesPerElement( mesh );
    GlobalOrdinal vertex_index;
    GlobalOrdinal vertex_ordinal;
    int destination_proc;
    Teuchos::Array<double> vertex_coords( d_vertex_dim );
    Teuchos::ArrayRCP<const double> mesh_coords =
	MeshTools<Mesh>::coordsView( mesh );
    Teuchos::ArrayRCP<const GlobalOrdinal> mesh_connectivity =
	MeshTools<Mesh>::connectivityView( mesh );
    for ( GlobalOrdinal n = 0; n < num_elements; ++n )
    {
	if ( elements_in_box[n] )
	{
	    for ( int i = 0; i < vertices_per_element; ++i )
	    {
		vertex_ordinal = mesh_connectivity[ i*num_elements + n ];
		vertex_index = vertex_indices.find( vertex_ordinal )->second;
		for ( int d = 0; d < d_vertex_dim; ++d )
		{
		    vertex_coords[d] = 
			mesh_coords[ d*num_vertices + vertex_index ];
		}
		destination_proc = d_rcb->getDestinationProc( vertex_coords );
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
    Teuchos::ArrayView<const int> from_images = 
	element_distributor.getImagesFrom();
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
    testInvariant( element_src_procs.size() == num_import_elements );
        
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

    // Now get the destination procs for all the vertices. This will be the
    // same destination procs as all of their parent elements. Therefore,
    // vertices may then also have to go to multiple procs because of this and
    // these procs may be different than their original RCB procs.
    GlobalOrdinal element_index;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
	export_elements_iterator;
    Teuchos::Array<int>::const_iterator export_element_procs_iterator;
    Teuchos::Array< std::set<int> > export_vertex_procs_set( num_vertices );
    for ( export_elements_iterator = export_elements.begin(),
     export_element_procs_iterator = export_element_procs.begin();
	  export_elements_iterator != export_elements.end();
	  ++export_elements_iterator, ++export_element_procs_iterator )
    {
	element_index = 
	    element_indices.find( *export_elements_iterator )->second;

	for ( int i = 0; i < vertices_per_element; ++i )
	{
	    vertex_ordinal = 
		mesh_connectivity[ i*num_elements + element_index ];
	    vertex_index = vertex_indices.find( vertex_ordinal )->second;

	    export_vertex_procs_set[ vertex_index ].insert( 
		*export_element_procs_iterator );
	}
    }
    export_elements.clear();    
    export_element_procs.clear();
    vertex_indices.clear();
    element_indices.clear();

    // Unroll the vector of sets into two vectors; one containing the vertex
    // ordinal and the other containing the corresponding vertex destination.
    Teuchos::Array<GlobalOrdinal> export_vertices;
    Teuchos::Array<int> export_vertex_procs;
    Teuchos::Array< std::set<int> >::const_iterator vertex_vec_iterator;
    std::set<int>::const_iterator vertex_proc_set_iterator;
    for ( vertex_vec_iterator = export_vertex_procs_set.begin(),
       export_vertex_iterator = MT::verticesBegin( mesh );
	  vertex_vec_iterator != export_vertex_procs_set.end();
	  ++vertex_vec_iterator, ++export_vertex_iterator )
    {
	for ( vertex_proc_set_iterator = vertex_vec_iterator->begin();
	      vertex_proc_set_iterator != vertex_vec_iterator->end();
	      ++vertex_proc_set_iterator )
	{
	    export_vertices.push_back( *export_vertex_iterator );
	    export_vertex_procs.push_back( *vertex_proc_set_iterator );
	}
    }
    export_vertex_procs_set.clear();

    // Now we know where the vertices need to go. Move the vertices to the
    // rendezvous decomposition through an inverse communciation operation.
    Tpetra::Distributor vertex_distributor( d_comm );
    Teuchos::ArrayView<int> export_vertex_procs_view = export_vertex_procs();
    GlobalOrdinal num_import_vertices = vertex_distributor.createFromSends(
	export_vertex_procs_view );
    Teuchos::ArrayView<const GlobalOrdinal> export_vertices_view = 
	export_vertices();
    Teuchos::Array<GlobalOrdinal> import_vertices( num_import_vertices );
    Teuchos::ArrayView<GlobalOrdinal> import_vertices_view = import_vertices();
    vertex_distributor.doPostsAndWaits( export_vertices_view, 1, 
					import_vertices_view );
    export_vertices.clear();
    export_vertex_procs.clear();

    // Next move these into the rendezvous vertex set so that we have a unique
    // list of the vertices.
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
	import_vertex_iterator;
    std::set<GlobalOrdinal> rendezvous_vertices_set;
    for ( import_vertex_iterator = import_vertices.begin();
	  import_vertex_iterator != import_vertices.end();
	  ++import_vertex_iterator )
    {
	rendezvous_vertices_set.insert( *import_vertex_iterator );
    }
    import_vertices.clear();

    // Finally put the vertices in a Teuchos::Array and get rid of the set.
    rendezvous_vertices.resize( rendezvous_vertices_set.size() );
    std::copy( rendezvous_vertices_set.begin(), rendezvous_vertices_set.end(),
	       rendezvous_vertices.begin() );
    rendezvous_vertices_set.clear();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_RENDEZVOUS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_Rendezvous_def.hpp
//---------------------------------------------------------------------------//
