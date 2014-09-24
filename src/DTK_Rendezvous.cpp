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
 * \file DTK_Rendezvous.cpp
 * \author Stuart R. Slattery
 * \brief Rendezvous definition.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <limits>

#include "DTK_Rendezvous.hpp"
#include "DTK_MeshTools.hpp"
#include "DTK_DBC.hpp"
#include "DTK_MeshTypes.hpp"
#include "DTK_PartitionerFactory.hpp"
#include "DTK_TopologyTools.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_as.hpp>

#include <Tpetra_Distributor.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_MultiVector.hpp>

#include <Intrepid_FieldContainer.hpp>

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
Rendezvous::Rendezvous( const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
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
Rendezvous::~Rendezvous()
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
void Rendezvous::build( const Teuchos::RCP<MeshManager>& mesh_manager,
			const CommIndexer& mesh_indexer,
			const BoundingBox& target_box )
{
    // Extract the mesh vertices and elements that are in the bounding
    // box. These are the pieces of the mesh that will be repartitioned.
    if ( Teuchos::nonnull(mesh_manager) ) 
    {
	getMeshInBox( mesh_manager );
    }

    // Construct the rendezvous partitioning for the mesh using the
    // vertices that are in the box.
    d_partitioner = PartitionerFactory::createMeshPartitioner( 
	d_comm, mesh_manager, d_dimension );
    DTK_ENSURE( Teuchos::nonnull(d_partitioner) );
    d_partitioner->partition( target_box );

    // Send the mesh in the box to the rendezvous decomposition and build the
    // mesh blocks.
    d_rendezvous_mesh_manager = sendMeshToRendezvous( mesh_manager, 
						      mesh_indexer );
    d_rendezvous_mesh_manager->buildIndexing();

    // Create a search tree in the rendezvous decomposition.
    d_element_tree = Teuchos::rcp( new ElementTree(d_rendezvous_mesh_manager) );
    DTK_ENSURE( Teuchos::nonnull(d_element_tree) );
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
Teuchos::Array<int> Rendezvous::procsContainingPoints(
    const Teuchos::ArrayRCP<double>& coords ) const
{
    Teuchos::Array<double> point( d_dimension );
    MeshId num_points = coords.size() / d_dimension;
    Teuchos::Array<int> destination_procs( num_points );
    for ( MeshId n = 0; n < num_points; ++n )
    {
	for ( int d = 0; d < d_dimension; ++d )
	{
	    point[d] = coords[ d*num_points + n ];
	}
	destination_procs[n] = d_partitioner->getPointDestinationProc( point() );
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
Teuchos::Array<Teuchos::Array<int> > Rendezvous::procsContainingBoxes( 
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
 * ordinal, std::numeric_limits<MeshId>::max(), for that point.
 *
 * \param element_src_procs The source procs that own the elements. Once proc
 * is provided for each element in the order that the elements were
 * provided. If a point is not found in an element, return an invalid element
 * source proc, -1, for that point.
 *
 * \param tolerance Absolute tolerance for point searching. Will be used when
 * checking the reference cell ( and is therefore absolute ).
 */
void Rendezvous::elementsContainingPoints( 
    const Teuchos::ArrayRCP<double>& coords,
    Teuchos::Array<MeshId>& elements,
    Teuchos::Array<int>& element_src_procs,
    double tolerance ) const
{
    Teuchos::Array<double> point( d_dimension );
    MeshId element_ordinal;
    MeshId num_points = coords.size() / d_dimension;
    bool found_point;
    elements.resize( num_points );
    element_src_procs.resize( num_points );
    for ( MeshId n = 0; n < num_points; ++n )
    {
	for ( int d = 0; d < d_dimension; ++d )
	{
	    point[d] = coords[ d*num_points + n ];
	}

	found_point = 
	    d_element_tree->findPoint( point(), element_ordinal, tolerance );

	if ( found_point )
	{
	    elements[n] = element_ordinal;
	    element_src_procs[n] = 
		d_element_src_procs_map.find( element_ordinal )->second;
	}
	else
	{
	    elements[n] = std::numeric_limits<MeshId>::max();
	    element_src_procs[n] = -1;
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the native elements in the rendezvous decomposition that are in
 * each geometric object in a list.
 *
 * \param geometry The geometric objects to check for element inclusion.
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
template<class Geometry> 
void Rendezvous::elementsInGeometry(
    const Teuchos::Array<Geometry>& geometry,
    Teuchos::Array<Teuchos::Array<MeshId> >& elements,
    const double tolerance, bool all_vertices_for_inclusion ) const
{
    elements.resize( geometry.size() );

    typename Teuchos::Array<Geometry>::const_iterator geometry_iterator;
    typename Teuchos::Array<Teuchos::Array<MeshId> >::iterator 
	element_iterator;
    int num_blocks = d_rendezvous_mesh_manager->getNumBlocks();
    Teuchos::RCP<MeshBlock> block;
    Teuchos::ArrayRCP<double> block_vertex_coords;
    Teuchos::ArrayRCP<MeshId> block_element_gids;
    Teuchos::ArrayRCP<MeshId> block_connectivity;
    int block_verts_per_elem = 0;
    int block_num_elements = 0;
    int block_num_vertices = 0;
    MeshId vertex_gid = 0;
    int vertex_lid = 0;

    // Mesh elements are the outer loop. That way we only have to extract the
    // coordinates once.
    for ( int b = 0; b < num_blocks; ++b )
    {
	block = d_rendezvous_mesh_manager->getBlock( b );
	block_num_elements = block->numElements();
	block_vertex_coords = block->vertexCoordinates();
	block_element_gids = block->elementIds();
	block_connectivity = block->connectivity();
	block_num_vertices = block->numVertices();
	block_verts_per_elem = block->verticesPerElement();
	Intrepid::FieldContainer<double> node_coords( 
	    1, block_verts_per_elem, d_rendezvous_mesh_manager->dim() );
	for ( int e = 0; e < block_num_elements; ++e )
	{
	    // Extract the element node coordinates.
	    for ( int n = 0; n < block_verts_per_elem; ++n )
	    {
		vertex_gid = block_connectivity[ n*block_num_elements + e ];
		vertex_lid = d_rendezvous_mesh_manager->getLocalVertexId( 
		    b, vertex_gid );

		for ( int d = 0; d < d_rendezvous_mesh_manager->dim(); ++d )
		{
		    node_coords( 0, n, d ) = 
			block_vertex_coords[ block_num_vertices*d + vertex_lid ];
		}
	    }

	    // Check for inclusion in the geometry.
	    for ( geometry_iterator = geometry.begin(), 
		   element_iterator = elements.begin();
		  geometry_iterator != geometry.end();
		  ++geometry_iterator, ++element_iterator )
	    {
		if ( TopologyTools::elementInGeometry(
			 *geometry_iterator, node_coords, 
			 tolerance, all_vertices_for_inclusion) )
		{
		    element_iterator->push_back( block_element_gids[e] );
		}
	    }
	}
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
Teuchos::Array<int> Rendezvous::elementSourceProcs(
    const Teuchos::Array<MeshId>& elements )
{
    Teuchos::Array<int> source_procs( elements.size() );
    Teuchos::Array<int>::iterator source_proc_iterator;
    typename Teuchos::Array<MeshId>::const_iterator element_iterator;
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
void Rendezvous::getMeshInBox( const Teuchos::RCP<MeshManager>& mesh_manager )
{
    // Expand the box by a typical mesh element length in all directions plus
    // some tolerance. Doing this catches a few corner cases.
    double tol = 1.0e-4;
    Teuchos::Tuple<double,6> box_bounds = d_global_box.getBounds();
    double box_volume = d_global_box.volume( d_dimension );
    MeshId global_num_elements = mesh_manager->globalNumElements();
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> >::iterator block_iterator;
    for ( block_iterator = mesh_blocks.begin();
	  block_iterator != mesh_blocks.end();
	  ++block_iterator )
    {
	// Setup.
	vertices_in_box.clear();
	elements_in_box.clear();
	int block_id = std::distance( mesh_blocks.begin(), 
				      block_iterator );

	// Create a map indexed by vertex global ordinal containing the actual
	// vertex ordinal location in the array and get all of the vertices
	// that are in the box.
	MeshId num_vertices = (*block_iterator)->numVertices();
	std::tr1::unordered_map<MeshId,MeshId> vertex_indices;
	Teuchos::ArrayRCP<MeshId> vertex_gids =
	    (*block_iterator)->vertexIds();
	Teuchos::Array<double> vertex_coords( d_dimension );
	Teuchos::ArrayRCP<double> mesh_coords =
	    (*block_iterator)->vertexCoordinates();
	for ( MeshId n = 0; n < num_vertices; ++n )
	{
	    // Set the mapping.
	    vertex_indices[ vertex_gids[n] ] = n;

	    // Get all of the vertices that are in the box. 
	    for ( int d = 0; d < d_dimension; ++d )
	    {
		vertex_coords[d] = mesh_coords[ d*num_vertices + n ];
	    }
	    vertices_in_box.push_back( 
		d_global_box.pointInBox( vertex_coords ) );
	}
	DTK_CHECK( Teuchos::as<MeshId>(vertices_in_box.size()) 
		   == num_vertices );

	// For those vertices that are in the box, get the elements that they
	// construct. These elements are in the box.
	MeshId num_elements = (*block_iterator)->numElements();
	int vertices_per_element = (*block_iterator)->verticesPerElement();
	MeshId vertex_index;
	MeshId vertex_ordinal;
	int this_element_in_box;
	Teuchos::ArrayRCP<MeshId> mesh_connectivity = 
	    (*block_iterator)->connectivity();
	for ( MeshId n = 0; n < num_elements; ++n )
	{
	    this_element_in_box = 0;
	    for ( int i = 0; i < vertices_per_element; ++i )
	    {
		vertex_ordinal = mesh_connectivity[ i*num_elements + n ];
		vertex_index = vertex_indices.find( vertex_ordinal )->second;
		if ( vertices_in_box[ vertex_index ] )
		{
		    this_element_in_box = 1;
		    break;
		}
	    }
	    elements_in_box.push_back( this_element_in_box );

	    // Get the vertices that belong to the elements in the box, but
	    // are not necessarily in the box themselves. These will also be
	    // used in partitioning.
	    if ( this_element_in_box )
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
	DTK_CHECK( Teuchos::as<MeshId>(elements_in_box.size())
		   == num_elements );

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
Teuchos::RCP<MeshManager> Rendezvous::sendMeshToRendezvous( 
    const Teuchos::RCP<MeshManager>& mesh_manager,
    const CommIndexer& mesh_indexer )
{
    // Set a value for mesh existence.
    bool mesh_exists = Teuchos::nonnull( mesh_manager );

    // Setup a mesh indexer and the mesh blocks.
    int num_mesh_blocks = 0;
    if ( mesh_exists )
    {
	num_mesh_blocks = mesh_manager->getNumBlocks();
    }
    Teuchos::broadcast<int,int>( *d_comm, mesh_indexer.l2g(0),
				 Teuchos::Ptr<int>(&num_mesh_blocks) );
    
    // Setup an array to store the mesh in the rendezvous decomposition.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> >
	rendezvous_block_containers( num_mesh_blocks );

    // Loop through the blocks and move them to the rendezvous decomposition.
    Teuchos::RCP<MeshBlock> current_block;
    for ( int block_id = 0; block_id < num_mesh_blocks; ++block_id )
    {
	// Get the block data.
	Teuchos::ArrayView<short int> active_block_elements;
	MeshId num_vertices = 0;
	Teuchos::ArrayRCP<MeshId> export_vertex_arcp(0,0);
	MeshId num_elements = 0;
	Teuchos::ArrayRCP<MeshId> export_element_arcp(0,0);
	Teuchos::ArrayRCP<double> export_coords_view(0,0);
	int max_verts_per_elem = 27;
	Teuchos::Array<int> block_packet(max_verts_per_elem+2,-1);
	Teuchos::ArrayRCP<MeshId> export_conn_view(0,0);
	int topo_int = 0;
	if ( mesh_exists )
	{
	    current_block = mesh_manager->getBlock( block_id );
	    active_block_elements = mesh_manager->getActiveElements( block_id );
	    num_vertices = current_block->numVertices();
	    export_vertex_arcp = current_block->vertexIds();
	    num_elements = current_block->numElements();
	    export_element_arcp = current_block->elementIds();
	    export_coords_view = current_block->vertexCoordinates();
	    export_conn_view = current_block->connectivity();
	    topo_int = static_cast<int>( current_block->elementTopology() );
	    block_packet[0] = current_block->verticesPerElement();
	    block_packet[1] = topo_int;
	    std::copy( current_block->permutation().begin(),
		       current_block->permutation().end(),
		       &block_packet[2] );
	}

	// Broadcast the number of vertices per element, the element
	// topology, and the permutation list.
	Teuchos::broadcast<int,int>( *d_comm, mesh_indexer.l2g(0), block_packet() );

	// Extract the number of vertices per element.
	int verts_per_elem = block_packet[0];

	// Extract the block topology.
	topo_int = block_packet[1];
	DTK_ElementTopology block_topology = 
	    static_cast<DTK_ElementTopology>( topo_int );

	// Extract the permutation list.
	DTK_CHECK( verts_per_elem <= max_verts_per_elem );
	Teuchos::ArrayRCP<int> permutation_list(verts_per_elem,0);
	std::copy( &block_packet[2], &block_packet[2] + verts_per_elem, 
		   permutation_list.begin() );
	block_packet.clear();

	// Setup the communication patterns for moving the mesh block to the
	// rendezvous decomposition. This will also move the vertex and element
	// global ordinals to the rendezvous decomposition.
	Teuchos::ArrayRCP<MeshId> rendezvous_vertices;
	Teuchos::ArrayRCP<MeshId> rendezvous_elements;
	setupImportCommunication( current_block, active_block_elements,
				  rendezvous_vertices, rendezvous_elements );

	// Setup export vertex map.
	Teuchos::ArrayView<const MeshId> export_vertex_view 
	    = export_vertex_arcp();
	Teuchos::RCP<const Tpetra::Map<int,MeshId> > export_vertex_map = 
	    Tpetra::createNonContigMap<int,MeshId>( 
		export_vertex_view, d_comm );
	DTK_CHECK( Teuchos::nonnull(export_vertex_map) );

	// Setup import vertex map.
	Teuchos::ArrayView<const MeshId> rendezvous_vertices_view = 
	    rendezvous_vertices();
	Teuchos::RCP<const Tpetra::Map<int,MeshId> > import_vertex_map = 
	    Tpetra::createNonContigMap<int,MeshId>(
		rendezvous_vertices_view, d_comm );
	DTK_CHECK( Teuchos::nonnull(import_vertex_map) );

	// Move the vertex coordinates to the rendezvous decomposition.
	Tpetra::Import<int,MeshId> vertex_importer( export_vertex_map, 
							   import_vertex_map );
	Teuchos::RCP< Tpetra::MultiVector<double,int,MeshId> > 
	    export_coords = Tpetra::createMultiVectorFromView( 
		export_vertex_map, export_coords_view, 
		num_vertices, d_dimension );
	Tpetra::MultiVector<double,int,MeshId> 
	    import_coords( import_vertex_map, d_dimension );
	import_coords.doImport( 
	    *export_coords, vertex_importer, Tpetra::INSERT );

	// Setup export element map.
	Teuchos::ArrayView<const MeshId> export_element_view = 
	    export_element_arcp();
	Teuchos::RCP<const Tpetra::Map<int,MeshId> > export_element_map = 
	    Tpetra::createNonContigMap<int,MeshId>(
		export_element_view, d_comm );
	DTK_CHECK( Teuchos::nonnull(export_element_map) );

	// Setup import element map.
	Teuchos::ArrayView<const MeshId> rendezvous_elements_view =
	    rendezvous_elements();
	Teuchos::RCP<const Tpetra::Map<int,MeshId> > import_element_map = 
	    Tpetra::createNonContigMap<int,MeshId>(
		rendezvous_elements_view, d_comm );
	DTK_CHECK( Teuchos::nonnull(import_element_map) );

	// Move the element connectivity to the rendezvous decomposition.
	Tpetra::Import<int,MeshId> element_importer( export_element_map, 
							    import_element_map );
	Teuchos::RCP<Tpetra::MultiVector<MeshId,int,MeshId> > 
	    export_conn = Tpetra::createMultiVectorFromView( 
		export_element_map, export_conn_view, 
		num_elements, verts_per_elem );
	Tpetra::MultiVector<MeshId,int,MeshId> import_conn( 
	    import_element_map, verts_per_elem );
	import_conn.doImport( *export_conn, element_importer, Tpetra::INSERT );

	// Construct the mesh block container from the collected data,
	// effectively wrapping it with mesh traits.
	rendezvous_block_containers[ block_id ] = 
	    Teuchos::rcp( new MeshContainer( 
			      d_dimension,
			      rendezvous_vertices,
			      import_coords.get1dViewNonConst(),
			      block_topology,
			      verts_per_elem,
			      rendezvous_elements,
			      import_conn.get1dViewNonConst(),
			      permutation_list ) );
    }

    // Build the rendezvous mesh manager from the rendezvous mesh blocks.
    return Teuchos::rcp( 
	new MeshManager( rendezvous_block_containers, d_comm, d_dimension) );
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
void Rendezvous::setupImportCommunication( 
    const Teuchos::RCP<MeshBlock>& mesh,
    const Teuchos::ArrayView<short int>& elements_in_box,
    Teuchos::ArrayRCP<MeshId>& rendezvous_vertices,
    Teuchos::ArrayRCP<MeshId>& rendezvous_elements )
{
    // Set a value for mesh existence.
    bool mesh_exists = Teuchos::nonnull( mesh );

    // Get the local block data and create a vertex index map and element
    // index map for access to connectivity data.
    MeshId num_vertices = 0;
    MeshId num_elements = 0;
    int vertices_per_element = 0;
    Teuchos::ArrayRCP<MeshId> mesh_vertices(0,0);
    Teuchos::ArrayRCP<MeshId> mesh_elements(0,0);
    Teuchos::ArrayRCP<MeshId> mesh_connectivity(0,0);
    Teuchos::ArrayRCP<MeshId>::const_iterator export_vertex_iterator;
    Teuchos::ArrayRCP<MeshId>::const_iterator export_element_iterator;
    std::tr1::unordered_map<MeshId,MeshId> vertex_indices;
    std::tr1::unordered_map<MeshId,MeshId> element_indices;
    MeshId array_index = 0;
    MeshId vertex_sum = 0;
    if ( mesh_exists )
    {
	// Get the block data.
	num_vertices = mesh->numVertices();
	num_elements = mesh->numElements();
	vertices_per_element = mesh->verticesPerElement();
	mesh_connectivity = mesh->connectivity();
	mesh_vertices = mesh->vertexIds();

	// Make the vertex map.
	for ( export_vertex_iterator = mesh_vertices.begin();
	      export_vertex_iterator != mesh_vertices.end();
	      ++export_vertex_iterator )
	{
	    vertex_indices[ *export_vertex_iterator ] = array_index;
	    ++array_index;
	}

	// Reset the index.
	array_index = 0;

	// Make the element map.
	for ( export_element_iterator = mesh_elements.begin();
	      export_element_iterator != mesh_elements.end();
	      ++export_element_iterator )
	{
	    element_indices[ *export_element_iterator ] = array_index;
	    ++array_index;
	}
    }

    // Get the destination procs for all local vertices in the global
    // bounding box.
    Teuchos::Array<int> vertex_procs =
	d_partitioner->getInputPointDestinationProcs( vertex_sum, num_vertices );
    DTK_CHECK( vertex_procs.size() == (int) num_vertices );
    vertex_sum += num_vertices;

    // Get destination procs for all local elements in the global bounding
    // box. The element will need to be sent to each partition that its
    // connecting vertices exist in. We'll make a unique destination proc set
    // for each element.
    Teuchos::Array< std::set<int> > export_element_procs_set( num_elements );
    MeshId vertex_index;
    MeshId vertex_ordinal;
    for ( MeshId n = 0; n < num_elements; ++n )
    {
	if ( elements_in_box[n] )
	{
	    for ( int i = 0; i < vertices_per_element; ++i )
	    {
		vertex_ordinal = mesh_connectivity[ i*num_elements + n ];
		vertex_index = vertex_indices.find( vertex_ordinal )->second;
		export_element_procs_set[n].insert(vertex_procs[vertex_index]);
	    }
	}
    }

    // Unroll the vector of sets into two vectors; one containing the element
    // ordinal and the other containing the corresponding element destination.
    Teuchos::Array<MeshId> export_elements;
    Teuchos::Array<int> export_element_procs;
    Teuchos::ArrayRCP<MeshId>::const_iterator element_iterator;
    Teuchos::Array<std::set<int> >::const_iterator element_vec_iterator;
    std::set<int>::const_iterator element_proc_set_iterator;
    if ( mesh_exists )
    {
	for ( element_vec_iterator = export_element_procs_set.begin(),
		  element_iterator = mesh_elements.begin();
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
    export_element_procs_set.clear();

    // Now we know where the elements need to go. Move the elements to the
    // rendezvous decomposition through an inverse communciation operation.
    Tpetra::Distributor element_distributor( d_comm );
    Teuchos::ArrayView<int> export_element_procs_view = export_element_procs();
    MeshId num_import_elements = element_distributor.createFromSends(
	export_element_procs_view );
    Teuchos::ArrayView<const MeshId> export_elements_view =
	export_elements();
    Teuchos::Array<MeshId> import_elements( num_import_elements );
    Teuchos::ArrayView<MeshId> import_elements_view = import_elements();
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
    DTK_CHECK( Teuchos::as<MeshId>(element_src_procs.size())
	       == num_import_elements );

    // Next, move these into the rendezvous element set so that we have a
    // unique list of the elements and build the rendezvous mesh element to
    // source proc map.
    std::set<MeshId> rendezvous_elements_set;
    for ( MeshId n = 0; n < num_import_elements; ++n )
    {
	if ( rendezvous_elements_set.insert( import_elements[n] ).second )
	{
	    d_element_src_procs_map[ import_elements[n] ] = 
		element_src_procs[n];
	}
    }
    import_elements.clear();
    element_src_procs.clear();

    // Finally put the elements in an array and get rid of the set.
    rendezvous_elements
	= Teuchos::ArrayRCP<MeshId>( rendezvous_elements_set.size() );
    std::copy( rendezvous_elements_set.begin(), rendezvous_elements_set.end(),
	       rendezvous_elements.begin() );
    rendezvous_elements_set.clear();

    // Now get the destination procs for all the vertices. This will be the
    // same destination procs as all of their parent elements. Therefore,
    // vertices may then also have to go to multiple procs because of this and
    // these procs may be different than their original rendezvous procs.
    MeshId element_index;
    typename Teuchos::Array<MeshId>::const_iterator 
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
    Teuchos::Array<MeshId> export_vertices;
    Teuchos::Array<int> export_vertex_procs;
    Teuchos::Array< std::set<int> >::const_iterator vertex_vec_iterator;
    std::set<int>::const_iterator vertex_proc_set_iterator;
    if ( mesh_exists )
    {
	for ( vertex_vec_iterator = export_vertex_procs_set.begin(),
	   export_vertex_iterator = mesh_vertices.begin();
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
    export_vertex_procs_set.clear();

    // Now we know where the vertices need to go. Move the vertices to the
    // rendezvous decomposition through an inverse communciation operation.
    Tpetra::Distributor vertex_distributor( d_comm );
    Teuchos::ArrayView<int> export_vertex_procs_view = export_vertex_procs();
    MeshId num_import_vertices = vertex_distributor.createFromSends(
	export_vertex_procs_view );
    Teuchos::ArrayView<const MeshId> export_vertices_view = 
	export_vertices();
    Teuchos::Array<MeshId> import_vertices( num_import_vertices );
    Teuchos::ArrayView<MeshId> import_vertices_view = import_vertices();
    vertex_distributor.doPostsAndWaits( export_vertices_view, 1, 
					import_vertices_view );
    export_vertices.clear();
    export_vertex_procs.clear();

    // Create a unique list of vertex ids.
    std::sort( import_vertices.begin(), import_vertices.end() );
    typename Teuchos::Array<MeshId>::iterator unique_vert_it =
	std::unique( import_vertices.begin(), import_vertices.end() );

    // Finally put the vertices in the output array.
    rendezvous_vertices = Teuchos::ArrayRCP<MeshId>( 
	std::distance(import_vertices.begin(), unique_vert_it) );
    std::copy( import_vertices.begin(), unique_vert_it,
	       rendezvous_vertices.begin() );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_Rendezvous.cpp
//---------------------------------------------------------------------------//
