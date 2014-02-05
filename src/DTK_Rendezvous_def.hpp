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
#include <limits>

#include "DTK_MeshTools.hpp"
#include "DTK_DBC.hpp"
#include "DTK_CommIndexer.hpp"
#include "DTK_MeshTypes.hpp"
#include "DTK_PartitionerFactory.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_as.hpp>

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
 * \param dimension The dimension of the rendezvous decomposition. We need
 * this here because the mesh we get may or may not exist on every
 * process. This prevents global communications.
 *
 * \param global_box The global bounding box inside of which the rendezvous
 * decomposition will be generated.
 */
template<class Mesh>
Rendezvous<Mesh>::Rendezvous( const RCP_Comm& comm,
			      const int dimension,
			      const BoundingBox& global_box )
    : d_comm( comm )
    , d_dimension( dimension )
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
 *
 * \param mesh_manager The mesh to repartition to the rendezvous
 * decomposition. A null argument is valid here as a mesh may not exist in its
 * original decomposition on every process in the global communicator. It
 * will, however, be redistributed across every process in the global
 * communicator. 
 */
template<class Mesh> 
void Rendezvous<Mesh>::build( const RCP_MeshManager& mesh_manager )
{
    // Extract the mesh vertices and elements that are in the bounding
    // box. These are the pieces of the mesh that will be repartitioned.
    if ( !mesh_manager.is_null() ) 
    {
	getMeshInBox( mesh_manager );
    }
    d_comm->barrier();

    // Construct the rendezvous partitioning for the mesh using the
    // vertices that are in the box.
    d_partitioner = PartitionerFactory::createMeshPartitioner( 
	d_comm, mesh_manager, d_dimension );
    DTK_ENSURE( Teuchos::nonnull(d_partitioner) );
    d_partitioner->partition();

    // Send the mesh in the box to the rendezvous decomposition and build the
    // mesh blocks.
    MeshManager<MeshContainerType> rendezvous_mesh_manager =
	sendMeshToRendezvous( mesh_manager );

    // Create a search tree in the rendezvous decomposition.
    d_search_tree = Teuchos::rcp(new ElementTree<Mesh>(rendezvous_mesh_manager));
    DTK_ENSURE( Teuchos::nonnull(d_search_tree) );
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
    Teuchos::Array<double> point( d_dimension );
    GlobalOrdinal num_points = coords.size() / d_dimension;
    Teuchos::Array<int> destination_procs( num_points );
    for ( GlobalOrdinal n = 0; n < num_points; ++n )
    {
	for ( int d = 0; d < d_dimension; ++d )
	{
	    point[d] = coords[ d*num_points + n ];
	}
	destination_procs[n] = d_partitioner->getPointDestinationProc( point );
    }

    return destination_procs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the rendezvous destination processes for a set of bounding
 * boxes. 
 *
 * \param boxes A list of boxes for which the rendezvous decomposition
 * destination procs are desired.
 *
 * \return A list of rendezvous decomposition procs for each box. A box may
 * have multiple procs that it spans.
 */
template<class Mesh>
Teuchos::Array<Teuchos::Array<int> > Rendezvous<Mesh>::procsContainingBoxes( 
    const Teuchos::Array<BoundingBox>& boxes ) const
{
    Teuchos::Array<Teuchos::Array<int> > box_procs( boxes.size() );
    Teuchos::Array<Teuchos::Array<int> >::iterator proc_iterator;
    Teuchos::Array<BoundingBox>::const_iterator box_iterator;
    for ( box_iterator = boxes.begin(), proc_iterator = box_procs.begin();
	  box_iterator != boxes.end();
	  ++box_iterator, ++proc_iterator )
    {
	*proc_iterator = d_partitioner->getBoxDestinationProcs( *box_iterator );
    }

    return box_procs;
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
 * ordinal, std::numeric_limits<GlobalOrdinal>::max(), for that point.
 *
 * \param element_src_procs The source procs that own the elements. Once proc
 * is provided for each element in the order that the elements were
 * provided. If a point is not found in an element, return an invalid element
 * source proc, -1, for that point.
 *
 * \param tolerance Absolute tolerance for point searching. Will be used when
 * checking the reference cell ( and is therefore absolute ).
 */
template<class Mesh>
void Rendezvous<Mesh>::elementsContainingPoints( 
    const Teuchos::ArrayRCP<double>& coords,
    Teuchos::Array<GlobalOrdinal>& elements,
    Teuchos::Array<int>& element_src_procs,
    double tolerance ) const
{
    Teuchos::Array<double> point( d_dimension );
    GlobalOrdinal element_ordinal;
    GlobalOrdinal num_points = coords.size() / d_dimension;
    bool found_point;
    elements.resize( num_points );
    element_src_procs.resize( num_points );
    for ( GlobalOrdinal n = 0; n < num_points; ++n )
    {
	for ( int d = 0; d < d_dimension; ++d )
	{
	    point[d] = coords[ d*num_points + n ];
	}

	found_point = d_kdtree->findPoint( point, element_ordinal, tolerance );

	if ( found_point )
	{
	    elements[n] = element_ordinal;
	    element_src_procs[n] = 
		d_element_src_procs_map.find( element_ordinal )->second;
	}
	else
	{
	    elements[n] = std::numeric_limits<GlobalOrdinal>::max();
	    element_src_procs[n] = -1;
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the native elements in the rendezvous decomposition that are in
 * each bounding box in a list.
 *
 * \param boxes The boxes to get the elements for.
 *
 * \param elements The elements found in each of the boxes.
 *
 * \param element_src_procs The source procs owning the elements found in the
 * boxes.
 */
template<class Mesh>
void Rendezvous<Mesh>::elementsInBoxes(
    const Teuchos::Array<BoundingBox>& boxes,
    Teuchos::Array<Teuchos::Array<GlobalOrdinal> >& elements ) const
{
    elements.resize( boxes.size() );

    Teuchos::Array<BoundingBox>::const_iterator box_iterator;
    typename Teuchos::Array<Teuchos::Array<GlobalOrdinal> >::iterator 
	element_iterator;
    for ( box_iterator = boxes.begin(), element_iterator = elements.begin();
	  box_iterator != boxes.end();
	  ++box_iterator, ++element_iterator )
    {
	*element_iterator = d_rendezvous_mesh->elementsInBox( *box_iterator );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the native elements in the rendezvous decomposition that are in
 * each geometric object in a list.
 *
 * \param geometry The geometric object to get the elements for.
 *
 * \param elements The elements found in each of the boxes.
 *
 * \param tolerance Tolerance used for element vertex-in-geometry
 * checks.
 *
 * \param all_vertices_for_inclusion Flag for element-in-geometry
 * inclusion. If set to true, all of an element's vertices are required to
 * reside within a geometry within the geometric tolerance in order to be
 * considered a member of that geometry's conformal mesh. If set to false,
 * only one of an element's vertices must be contained within the geometric
 * tolerance of the geometry in order to be considered a member of that
 * geometry's conformal mesh.
 */
template<class Mesh>
template<class Geometry> 
void Rendezvous<Mesh>::elementsInGeometry(
    const Teuchos::Array<Geometry>& geometry,
    Teuchos::Array<Teuchos::Array<GlobalOrdinal> >& elements,
    const double tolerance, bool all_vertices_for_inclusion ) const
{
    elements.resize( geometry.size() );

    typename Teuchos::Array<Geometry>::const_iterator geometry_iterator;
    typename Teuchos::Array<Teuchos::Array<GlobalOrdinal> >::iterator 
	element_iterator;
    for ( geometry_iterator = geometry.begin(), 
	   element_iterator = elements.begin();
	  geometry_iterator != geometry.end();
	  ++geometry_iterator, ++element_iterator )
    {
	*element_iterator = 
	    d_rendezvous_mesh->elementsInGeometry( *geometry_iterator,
						   tolerance,
						   all_vertices_for_inclusion );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief For a list of elements in the rendezvous decomposition, the their
 * source procs.
 *
 * \param elements The elements to get source procs for.
 *
 * \return The source procs for the elements.
 */
template<class Mesh>
Teuchos::Array<int> Rendezvous<Mesh>::elementSourceProcs(
    const Teuchos::Array<GlobalOrdinal>& elements )
{
    Teuchos::Array<int> source_procs( elements.size() );
    Teuchos::Array<int>::iterator source_proc_iterator;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator element_iterator;
    for ( element_iterator = elements.begin(),
      source_proc_iterator = source_procs.begin();
	  element_iterator != elements.end();
	  ++element_iterator, ++source_proc_iterator )
    {
	*source_proc_iterator = 
	    d_element_src_procs_map.find( *element_iterator )->second;
    }

    return source_procs;
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
    double box_volume = d_global_box.volume( d_dimension );
    GlobalOrdinal global_num_elements = mesh_manager->globalNumElements();
    double pow = 1 / d_dimension;
    double typical_length = std::pow( box_volume / global_num_elements, pow );
    for ( int d = 0; d < d_dimension; ++d )
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
	    MeshTools<Mesh>::numVertices( *(*block_iterator) );
	std::tr1::unordered_map<GlobalOrdinal,GlobalOrdinal> vertex_indices;
	typename MT::const_vertex_iterator vertex_iterator;
	GlobalOrdinal array_index = 0;
	for ( vertex_iterator = MT::verticesBegin( *(*block_iterator) );
	      vertex_iterator != MT::verticesEnd( *(*block_iterator) );
	      ++vertex_iterator )
	{
	    vertex_indices[ *vertex_iterator ] = array_index;
	    ++array_index;
	}

	// Get all of the vertices that are in the box. 
	Teuchos::Array<double> vertex_coords( d_dimension );
	Teuchos::ArrayRCP<const double> mesh_coords =
	    MeshTools<Mesh>::coordsView( *(*block_iterator) );
	for ( GlobalOrdinal n = 0; n < num_vertices; ++n )
	{
	    for ( int d = 0; d < d_dimension; ++d )
	    {
		vertex_coords[d] = mesh_coords[ d*num_vertices + n ];
	    }
	    vertices_in_box.push_back( 
		d_global_box.pointInBox( vertex_coords ) );
	}
	DTK_CHECK( Teuchos::as<GlobalOrdinal>(vertices_in_box.size()) 
		       == num_vertices );

	// For those vertices that are in the box, get the elements that they
	// construct. These elements are in the box.
	GlobalOrdinal num_elements = 
	    MeshTools<Mesh>::numElements( *(*block_iterator) );
	int vertices_per_element = MT::verticesPerElement( *(*block_iterator) );
	GlobalOrdinal vertex_index;
	GlobalOrdinal vertex_ordinal;
	int this_element_in_box;
	Teuchos::ArrayRCP<const GlobalOrdinal> mesh_connectivity = 
	    MeshTools<Mesh>::connectivityView( *(*block_iterator) );
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
	DTK_CHECK( Teuchos::as<GlobalOrdinal>(elements_in_box.size())
		       == num_elements );

	// Get the vertices that belong to the elements in the box, but are
	// not necessarily in the box themselves. These will also be used in
	// partitioning.
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
    // Set a value for mesh existence.
    bool mesh_exists = true;
    if ( mesh_manager.is_null() ) mesh_exists = false;

    // Setup a mesh indexer.
    RCP_Comm mesh_comm;
    if ( mesh_exists )
    {
	mesh_comm = mesh_manager->comm();
    }
    d_comm->barrier();
    CommIndexer mesh_indexer( d_comm, mesh_comm );

    // Setup the mesh blocks.
    int num_mesh_blocks;
    if ( mesh_exists )
    {
	num_mesh_blocks = mesh_manager->getNumBlocks();
    }
    d_comm->barrier();
    Teuchos::broadcast<int,int>( *d_comm, mesh_indexer.l2g(0),
				 Teuchos::Ptr<int>(&num_mesh_blocks) );
    
    // Setup an array to store the mesh in the rendezvous decomposition.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshContainerType> >
	rendezvous_block_containers( num_mesh_blocks );

    // Loop through the blocks and move them to the rendezvous decomposition.
    Teuchos::RCP<Mesh> current_block;
    for ( int block_id = 0; block_id < num_mesh_blocks; ++block_id )
    {
	// Get the current mesh block we are operating on.
	if ( mesh_exists )
	{
	    current_block = mesh_manager->getBlock( block_id );
	}
	d_comm->barrier();

	// Setup the communication patterns for moving the mesh block to the
	// rendezvous decomposition. This will also move the vertex and element
	// global ordinals to the rendezvous decomposition.
	Teuchos::Array<GlobalOrdinal> rendezvous_vertices;
	Teuchos::Array<GlobalOrdinal> rendezvous_elements;
	Teuchos::ArrayView<short int> active_block_elements;
	if ( mesh_exists )
	{
	    active_block_elements = mesh_manager->getActiveElements( block_id );
	}
	d_comm->barrier();
	setupImportCommunication( current_block, active_block_elements,
				  rendezvous_vertices, rendezvous_elements );

	// Setup export vertex map.
	GlobalOrdinal num_vertices = 0;
	Teuchos::ArrayRCP<GlobalOrdinal> export_vertex_arcp(0,0);
	if ( mesh_exists )
	{
	    num_vertices = MeshTools<Mesh>::numVertices( *current_block );
	    export_vertex_arcp = 
		MeshTools<Mesh>::verticesNonConstView( *current_block );
	}
	d_comm->barrier();
	Teuchos::ArrayView<const GlobalOrdinal> export_vertex_view 
	    = export_vertex_arcp();
	RCP_TpetraMap export_vertex_map = 
	    Tpetra::createNonContigMap<int,GlobalOrdinal>( 
		export_vertex_view, d_comm );
	DTK_CHECK( !export_vertex_map.is_null() );

	// Setup import vertex map.
	Teuchos::ArrayView<const GlobalOrdinal> rendezvous_vertices_view = 
	    rendezvous_vertices();
	RCP_TpetraMap import_vertex_map = 
	    Tpetra::createNonContigMap<int,GlobalOrdinal>(
		rendezvous_vertices_view, d_comm );
	DTK_CHECK( !import_vertex_map.is_null() );

	// Setup export element map.
	GlobalOrdinal num_elements = 0;
	Teuchos::ArrayRCP<GlobalOrdinal> export_element_arcp(0,0);
	if ( mesh_exists )
	{
	    num_elements = MeshTools<Mesh>::numElements( *current_block );
	    export_element_arcp =
		MeshTools<Mesh>::elementsNonConstView( *current_block );
	}
	d_comm->barrier();
	Teuchos::ArrayView<const GlobalOrdinal> export_element_view = 
	    export_element_arcp();
	RCP_TpetraMap export_element_map = 
	    Tpetra::createNonContigMap<int,GlobalOrdinal>(
		export_element_view, d_comm );
	DTK_CHECK( !export_element_map.is_null() );

	// Setup import element map.
	Teuchos::ArrayView<const GlobalOrdinal> rendezvous_elements_view =
	    rendezvous_elements();
	RCP_TpetraMap import_element_map = 
	    Tpetra::createNonContigMap<int,GlobalOrdinal>(
		rendezvous_elements_view, d_comm );
	DTK_CHECK( !import_element_map.is_null() );

	// Setup importers.
	Tpetra::Import<int,GlobalOrdinal> vertex_importer( export_vertex_map, 
							   import_vertex_map );
	Tpetra::Import<int,GlobalOrdinal> element_importer( export_element_map, 
							    import_element_map );

	// Move the vertex coordinates to the rendezvous decomposition.
	Teuchos::ArrayRCP<double> export_coords_view(0,0);
	if ( mesh_exists )
	{
	    export_coords_view =
		MeshTools<Mesh>::coordsNonConstView( *current_block );
	}
	d_comm->barrier();
	Teuchos::RCP< Tpetra::MultiVector<double,int,GlobalOrdinal> > 
	    export_coords = Tpetra::createMultiVectorFromView( 
		export_vertex_map, export_coords_view, 
		num_vertices, d_dimension );
	Tpetra::MultiVector<double,int,GlobalOrdinal> 
	    import_coords( import_vertex_map, d_dimension );
	import_coords.doImport( 
	    *export_coords, vertex_importer, Tpetra::INSERT );

	// Move the element connectivity to the rendezvous decomposition.
	int vertices_per_element;
	if ( mesh_exists )
	{
	    vertices_per_element = MT::verticesPerElement( *current_block );
	}
	d_comm->barrier();
	Teuchos::broadcast<int,int>( *d_comm, mesh_indexer.l2g(0),
				     Teuchos::Ptr<int>(&vertices_per_element) );

	Teuchos::ArrayRCP<GlobalOrdinal> export_conn_view(0,0);
	if ( mesh_exists ) 
	{
	    export_conn_view =
		MeshTools<Mesh>::connectivityNonConstView( *current_block );
	}
	d_comm->barrier();
	Teuchos::RCP<Tpetra::MultiVector<GlobalOrdinal,int,GlobalOrdinal> > 
	    export_conn = Tpetra::createMultiVectorFromView( 
		export_element_map, export_conn_view, 
		num_elements, vertices_per_element );
	Tpetra::MultiVector<GlobalOrdinal,int,GlobalOrdinal> import_conn( 
	    import_element_map, vertices_per_element );
	import_conn.doImport( *export_conn, element_importer, Tpetra::INSERT );

	// Broadcast the permutation list.
	Teuchos::ArrayRCP<int> permutation_list(vertices_per_element,0);
	if ( mesh_exists )
	{
	    permutation_list = 
		MeshTools<Mesh>::permutationNonConstView( *current_block );
	}
	d_comm->barrier();
	Teuchos::broadcast<int,int>( *d_comm, mesh_indexer.l2g(0),
				     vertices_per_element, 
				     &permutation_list[0] );

	// broadcast the element topology.
	int block_topology = 0;
	if ( mesh_exists )
	{
	    block_topology = 
		static_cast<int>(MT::elementTopology( *current_block ));
	}
	d_comm->barrier();
	Teuchos::broadcast<int,int>( *d_comm, mesh_indexer.l2g(0),
				     Teuchos::Ptr<int>(&block_topology) );

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

	rendezvous_block_containers[ block_id ] = 
	    Teuchos::rcp( new MeshContainerType( 
			      d_dimension,
			      rendezvous_vertices_array,
			      import_coords.get1dView(),
			      static_cast<DTK_ElementTopology>(block_topology),
			      vertices_per_element,
			      rendezvous_elements_array,
			      import_conn.get1dView(),
			      permutation_list ) );
    }

    // Build the rendezvous mesh manager from the rendezvous mesh blocks.
    return MeshManager<MeshContainerType>( rendezvous_block_containers,
					   d_comm,
					   d_dimension );
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
    const Teuchos::RCP<Mesh>& mesh,
    const Teuchos::ArrayView<short int>& elements_in_box,
    Teuchos::Array<GlobalOrdinal>& rendezvous_vertices,
    Teuchos::Array<GlobalOrdinal>& rendezvous_elements )
{
    // Set a value for mesh existence.
    bool mesh_exists = true;
    if ( mesh.is_null() ) mesh_exists = false;

    // Create a vertex index map for logarithmic time access to connectivity
    // data. 
    typename MT::const_vertex_iterator export_vertex_iterator;
    std::tr1::unordered_map<GlobalOrdinal,GlobalOrdinal> vertex_indices;
    GlobalOrdinal array_index = 0;
    if ( mesh_exists )
    {
	for ( export_vertex_iterator = MT::verticesBegin( *mesh );
	      export_vertex_iterator != MT::verticesEnd( *mesh );
	      ++export_vertex_iterator )
	{
	    vertex_indices[ *export_vertex_iterator ] = array_index;
	    ++array_index;
	}
    }
    d_comm->barrier();

    // Create a element index map for logarithmic time access to connectivity
    // data. 
    typename MT::const_element_iterator export_element_iterator;
    std::tr1::unordered_map<GlobalOrdinal,GlobalOrdinal> element_indices;
    array_index = 0;
    if ( mesh_exists )
    {
	for ( export_element_iterator = MT::elementsBegin( *mesh );
	      export_element_iterator != MT::elementsEnd( *mesh );
	      ++export_element_iterator )
	{
	    element_indices[ *export_element_iterator ] = array_index;
	    ++array_index;
	}
    }
    d_comm->barrier();

    // Get destination procs for all local elements in the global bounding
    // box. The element will need to be sent to each partition that its
    // connecting vertices exist in. We'll make a unique destination proc set
    // for each element.
    GlobalOrdinal num_vertices = 0;
    GlobalOrdinal num_elements = 0;
    int vertices_per_element = 0;
    if ( mesh_exists )
    {
	num_vertices = MeshTools<Mesh>::numVertices( *mesh );
	num_elements = MeshTools<Mesh>::numElements( *mesh );
	vertices_per_element = MT::verticesPerElement( *mesh );
    }
    d_comm->barrier();

    Teuchos::Array< std::set<int> > export_element_procs_set( num_elements );
    GlobalOrdinal vertex_index;
    GlobalOrdinal vertex_ordinal;
    int destination_proc;
    Teuchos::Array<double> vertex_coords( d_dimension );

    Teuchos::ArrayRCP<double> mesh_coords(0,0);
    Teuchos::ArrayRCP<GlobalOrdinal> mesh_connectivity(0,0);
    if ( mesh_exists )
    {
	mesh_coords = MeshTools<Mesh>::coordsNonConstView( *mesh );
	mesh_connectivity = MeshTools<Mesh>::connectivityNonConstView( *mesh );
    }
    d_comm->barrier();

    for ( GlobalOrdinal n = 0; n < num_elements; ++n )
    {
	if ( elements_in_box[n] )
	{
	    for ( int i = 0; i < vertices_per_element; ++i )
	    {
		vertex_ordinal = mesh_connectivity[ i*num_elements + n ];
		vertex_index = vertex_indices.find( vertex_ordinal )->second;
		for ( int d = 0; d < d_dimension; ++d )
		{
		    vertex_coords[d] = 
			mesh_coords[ d*num_vertices + vertex_index ];
		}
		destination_proc = 
		    d_partitioner->getPointDestinationProc( vertex_coords );
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
    if ( mesh_exists )
    {
	for ( element_vec_iterator = export_element_procs_set.begin(),
		  element_iterator = MT::elementsBegin( *mesh );
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
    }
    d_comm->barrier();

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
    DTK_CHECK( Teuchos::as<GlobalOrdinal>(element_src_procs.size())
		   == num_import_elements );
        
    // Next, move these into the rendezvous element set so that we have a
    // unique list of the elements and build the rendezvous mesh element to
    // source proc map.
    std::set<GlobalOrdinal> rendezvous_elements_set;
    for ( GlobalOrdinal n = 0; n < num_import_elements; ++n )
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
    // these procs may be different than their original rendezvous procs.
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
    if ( mesh_exists )
    {
	for ( vertex_vec_iterator = export_vertex_procs_set.begin(),
	   export_vertex_iterator = MT::verticesBegin( *mesh );
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
    }
    d_comm->barrier();

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
