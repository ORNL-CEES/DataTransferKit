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
 * \file DTK_GeometryRendezvous_def.hpp
 * \author Stuart R. Slattery
 * \brief GeometryRendezvous definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_GEOMETRYRENDEZVOUS_DEF_HPP
#define DTK_GEOMETRYRENDEZVOUS_DEF_HPP

#include <set>
#include <algorithm>
#include <limits>

#include "DTK_Assertion.hpp"
#include "DTK_CommIndexer.hpp"
#include "DTK_PartitionerFactory.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_as.hpp>

#include <Tpetra_Distributor.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Tpetra_MultiVector_def.hpp>

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
 * this here because the geometry we get may or may not exist on every
 * process. This prevents global communications.
 *
 * \param global_box The global bounding box inside of which the rendezvous
 * decomposition will be generated.
 */
template<class Geometry, class GlobalOrdinal>
GeometryRendezvous<Geometry,GlobalOrdinal>::GeometryRendezvous( 
    const RCP_Comm& comm,
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
template<class Geometry, class GlobalOrdinal>
GeometryRendezvous<Geometry,GlobalOrdinal>::~GeometryRendezvous()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Build the rendezvous decomposition.
 *
 * \param geometry_manager The geometry to repartition to the rendezvous
 * decomposition. A null argument is valid here as a geometry may not exist in its
 * original decomposition on every process in the global communicator. It
 * will, however, be redistributed across every process in the global
 * communicator. 
 */
template<class Geometry, class GlobalOrdinal> 
void GeometryRendezvous<Geometry,GlobalOrdinal>::build( 
    const RCP_GeometryManager& geometry_manager )
{
    // Extract the geometry vertices and elements that are in the bounding
    // box. These are the pieces of the geometry that will be repartitioned.
    if ( !geometry_manager.is_null() ) 
    {
	getGeometryInBox( geometry_manager );
    }

    // Construct the rendezvous partitioning for the geometry using the
    // vertices that are in the box.
    d_partitioner = PartitionerFactory::createGeometryPartitioner( 
	d_comm, geometry_manager, d_dimension );
    testPostcondition( !d_partitioner.is_null() );
    d_partitioner->partition();

    // Send the geometry in the box to the rendezvous decomposition.
    sendGeometryToRendezvous( geometry_manager );
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
template<class Geometry, class GlobalOrdinal>
Teuchos::Array<int> 
GeometryRendezvous<Geometry,GlobalOrdinal>::procsContainingPoints(
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
template<class Geometry, class GlobalOrdinal>
Teuchos::Array<Teuchos::Array<int> > 
GeometryRendezvous<Geometry,GlobalOrdinal>::procsContainingBoxes( 
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
 * \brief Get the native geometry ids in the rendezvous decomposition
 * containing a blocked list of coordinates also in the rendezvous
 * decomposition.
 * 
 * \param coords A blocked list of coordinates to search the geometry with. 

 * \param geometry An array of the geometries the points were found in. An
 * element will be returned for each point in the order they were provided
 * in. If a point is not found in an geometry, return an invalid geometry
 * ordinal, std::numeric_limits<GlobalOrdinal>::max(), for that point.
 *
 * \param geometry_src_procs The source procs that own the geometries. Once
 * proc is provided for each geometry in the order that the geometrys were
 * provided. If a point is not found in an geometry, return an invalid
 * geometry source proc, -1, for that point.
 */
template<class Geometry, class GlobalOrdinal>
void GeometryRendezvous<Geometry,GlobalOrdinal>::geometrysContainingPoints( 
    const Teuchos::ArrayRCP<double>& coords,
    Teuchos::Array<GlobalOrdinal>& geometry,
    Teuchos::Array<int>& geometry_src_procs ) const
{
    Teuchos::Array<double> point( d_dimension );
    GlobalOrdinal geometry_ordinal;
    GlobalOrdinal num_points = coords.size() / d_dimension;
    bool found_point;
    geometry.resize( num_points );
    geometry_src_procs.resize( num_points );
    for ( GlobalOrdinal n = 0; n < num_points; ++n )
    {
	for ( int d = 0; d < d_dimension; ++d )
	{
	    point[d] = coords[ d*num_points + n ];
	}

	found_point = d_kdtree->findPoint( point, geometry_ordinal );

	if ( found_point )
	{
	    geometry[n] = geometry_ordinal;
	    geometry_src_procs[n] = 
		d_geometry_src_procs_map.find( geometry_ordinal )->second;
	}
	else
	{
	    geometry[n] = std::numeric_limits<GlobalOrdinal>::max();
	    geometry_src_procs[n] = -1;
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief For a list of geometry in the rendezvous decomposition, the their
 * source procs.
 *
 * \param geometry The geometry to get source procs for.
 *
 * \return The source procs for the geometry.
 */
template<class Geometry, class GlobalOrdinal>
Teuchos::Array<int> 
GeometryRendezvous<Geometry,GlobalOrdinal>::geometrySourceProcs(
    const Teuchos::Array<GlobalOrdinal>& geometry )
{
    Teuchos::Array<int> source_procs( geometry.size() );
    Teuchos::Array<int>::iterator source_proc_iterator;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator geometry_iterator;
    for ( geometry_iterator = geometry.begin(),
      source_proc_iterator = source_procs.begin();
	  geometry_iterator != geometry.end();
	  ++geometry_iterator, ++source_proc_iterator )
    {
	*source_proc_iterator = 
	    d_geometry_src_procs_map.find( *geometry_iterator )->second;
    }

    return source_procs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Extract the geometry vertices and elements that are in a bounding box.
 *
 * \param geometry_manager The geometry to search the box with.
 */
template<class Geometry, class GlobalOrdinal>
GeometryRendezvous<Geometry,GlobalOrdinal>::getGeometryInBox( 
    const RCP_GeometryManager& geometry_manager )
{

}

//---------------------------------------------------------------------------//
/*!
 * \brief Send the geometry to the rendezvous decomposition.
 *
 * \param geometry_manager The geometry to send to the rendezvous decomposition.
 */
template<class Geometry, class GlobalOrdinal>
void GeometryRendezvous<Geometry,GlobalOrdinal>::sendGeometryToRendezvous( 
    const RCP_GeometryManager& geometry_manager )
{
    // Set a value for geometry existence.
    bool geometry_exists = true;
    if ( geometry_manager.is_null() ) geometry_exists = false;

    // Setup a geometry indexer.
    RCP_Comm geometry_comm;
    if ( geometry_exists )
    {
	geometry_comm = geometry_manager->comm();
    }
    d_comm->barrier();
    CommIndexer geometry_indexer( d_comm, geometry_comm );

    // Setup the geometry blocks.
    int num_geometry_blocks;
    if ( geometry_exists )
    {
	num_geometry_blocks = geometry_manager->getNumBlocks();
    }
    d_comm->barrier();
    Teuchos::broadcast<int,int>( *d_comm, geometry_indexer.l2g(0),
				 Teuchos::Ptr<int>(&num_geometry_blocks) );
    
    // Setup an array to store the geometry in the rendezvous decomposition.
    Teuchos::ArrayRCP<Teuchos::RCP<GeometryContainerType> >
	rendezvous_block_containers( num_geometry_blocks );

    // Loop through the blocks and move them to the rendezvous decomposition.
    Teuchos::RCP<Geometry> current_block;
    for ( int block_id = 0; block_id < num_geometry_blocks; ++block_id )
    {
	// Get the current geometry block we are operating on.
	if ( geometry_exists )
	{
	    current_block = geometry_manager->getBlock( block_id );
	}
	d_comm->barrier();

	// Setup the communication patterns for moving the geometry block to the
	// rendezvous decomposition. This will also move the vertex and element
	// global ordinals to the rendezvous decomposition.
	Teuchos::Array<GlobalOrdinal> rendezvous_vertices;
	Teuchos::Array<GlobalOrdinal> rendezvous_elements;
	Teuchos::ArrayView<short int> active_block_elements;
	if ( geometry_exists )
	{
	    active_block_elements = geometry_manager->getActiveElements( block_id );
	}
	d_comm->barrier();
	setupImportCommunication( current_block, active_block_elements,
				  rendezvous_vertices, rendezvous_elements );

	// Setup export vertex map.
	GlobalOrdinal num_vertices = 0;
	Teuchos::ArrayRCP<GlobalOrdinal> export_vertex_arcp(0,0);
	if ( geometry_exists )
	{
	    num_vertices = GeometryTools<Geometry>::numVertices( *current_block );
	    export_vertex_arcp = 
		GeometryTools<Geometry>::verticesNonConstView( *current_block );
	}
	d_comm->barrier();
	Teuchos::ArrayView<const GlobalOrdinal> export_vertex_view 
	    = export_vertex_arcp();
	RCP_TpetraMap export_vertex_map = 
	    Tpetra::createNonContigMap<GlobalOrdinal>( 
		export_vertex_view, d_comm );
	testInvariant( !export_vertex_map.is_null() );

	// Setup import vertex map.
	Teuchos::ArrayView<const GlobalOrdinal> rendezvous_vertices_view = 
	    rendezvous_vertices();
	RCP_TpetraMap import_vertex_map = 
	    Tpetra::createNonContigMap<GlobalOrdinal>(
		rendezvous_vertices_view, d_comm );
	testInvariant( !import_vertex_map.is_null() );

	// Setup export element map.
	GlobalOrdinal num_elements = 0;
	Teuchos::ArrayRCP<GlobalOrdinal> export_element_arcp(0,0);
	if ( geometry_exists )
	{
	    num_elements = GeometryTools<Geometry>::numElements( *current_block );
	    export_element_arcp =
		GeometryTools<Geometry>::elementsNonConstView( *current_block );
	}
	d_comm->barrier();
	Teuchos::ArrayView<const GlobalOrdinal> export_element_view = 
	    export_element_arcp();
	RCP_TpetraMap export_element_map = 
	    Tpetra::createNonContigMap<GlobalOrdinal>(
		export_element_view, d_comm );
	testInvariant( !export_element_map.is_null() );

	// Setup import element map.
	Teuchos::ArrayView<const GlobalOrdinal> rendezvous_elements_view =
	    rendezvous_elements();
	RCP_TpetraMap import_element_map = 
	    Tpetra::createNonContigMap<GlobalOrdinal>(
		rendezvous_elements_view, d_comm );
	testInvariant( !import_element_map.is_null() );

	// Setup importers.
	Tpetra::Import<GlobalOrdinal> vertex_importer( export_vertex_map, 
						       import_vertex_map );
	Tpetra::Import<GlobalOrdinal> element_importer( export_element_map, 
							import_element_map );

	// Move the vertex coordinates to the rendezvous decomposition.
	Teuchos::ArrayRCP<double> export_coords_view(0,0);
	if ( geometry_exists )
	{
	    export_coords_view =
		GeometryTools<Geometry>::coordsNonConstView( *current_block );
	}
	d_comm->barrier();
	Teuchos::RCP< Tpetra::MultiVector<double,GlobalOrdinal> > 
	    export_coords = Tpetra::createMultiVectorFromView( 
		export_vertex_map, export_coords_view, 
		num_vertices, d_dimension );
	Tpetra::MultiVector<double,GlobalOrdinal> 
	    import_coords( import_vertex_map, d_dimension );
	import_coords.doImport( 
	    *export_coords, vertex_importer, Tpetra::INSERT );

	// Move the element connectivity to the rendezvous decomposition.
	int vertices_per_element;
	if ( geometry_exists )
	{
	    vertices_per_element = MT::verticesPerElement( *current_block );
	}
	d_comm->barrier();
	Teuchos::broadcast<int,int>( *d_comm, geometry_indexer.l2g(0),
				     Teuchos::Ptr<int>(&vertices_per_element) );

	Teuchos::ArrayRCP<GlobalOrdinal> export_conn_view(0,0);
	if ( geometry_exists ) 
	{
	    export_conn_view =
		GeometryTools<Geometry>::connectivityNonConstView( *current_block );
	}
	d_comm->barrier();
	Teuchos::RCP<Tpetra::MultiVector<GlobalOrdinal,GlobalOrdinal> > 
	    export_conn = Tpetra::createMultiVectorFromView( 
		export_element_map, export_conn_view, 
		num_elements, vertices_per_element );
	Tpetra::MultiVector<GlobalOrdinal,GlobalOrdinal> import_conn( 
	    import_element_map, vertices_per_element );
	import_conn.doImport( *export_conn, element_importer, Tpetra::INSERT );

	// Broadcast the permutation list.
	Teuchos::ArrayRCP<int> permutation_list(vertices_per_element,0);
	if ( geometry_exists )
	{
	    permutation_list = 
		GeometryTools<Geometry>::permutationNonConstView( *current_block );
	}
	d_comm->barrier();
	Teuchos::broadcast<int,int>( *d_comm, geometry_indexer.l2g(0),
				     vertices_per_element, 
				     &permutation_list[0] );

	// broadcast the element topology.
	int block_topology = 0;
	if ( geometry_exists )
	{
	    block_topology = 
		static_cast<int>(MT::elementTopology( *current_block ));
	}
	d_comm->barrier();
	Teuchos::broadcast<int,int>( *d_comm, geometry_indexer.l2g(0),
				     Teuchos::Ptr<int>(&block_topology) );

	// Construct the geometry block container from the collected data,
	// effectively wrapping it with geometry traits.
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
	    Teuchos::rcp( new GeometryContainerType( 
			      d_dimension,
			      rendezvous_vertices_array,
			      import_coords.get1dView(),
			      static_cast<DTK_ElementTopology>(block_topology),
			      vertices_per_element,
			      rendezvous_elements_array,
			      import_conn.get1dView(),
			      permutation_list ) );
    }

    // Build the rendezvous geometry manager from the rendezvous geometry blocks.
    return GeometryManagerType( 
	rendezvous_block_containers, d_comm, d_dimension );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_GEOMETRYRENDEZVOUS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_GeometryRendezvous_def.hpp
//---------------------------------------------------------------------------//
