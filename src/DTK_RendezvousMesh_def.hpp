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
 * \file DTK_RendezvousMesh_def.hpp
 * \author Stuart R. Slattery
 * \brief Concrete Moab mesh template definitions.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_RENDEZVOUSMESH_DEF_HPP
#define DTK_RENDEZVOUSMESH_DEF_HPP

#include <algorithm>

#include "DTK_MeshTools.hpp"
#include "DTK_BoundingBox.hpp"
#include "DTK_TopologyTools.hpp"
#include "DTK_DBC.hpp"
#include "DataTransferKit_config.hpp"

#include <MBCore.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_as.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 * 
 * \param moab The Moab interface to build the RendezvousMesh with.
 *
 * \param ordinal_map A map relating the client element global ordinals to the
 * Moab element handles.
 */
template<typename GlobalOrdinal>
RendezvousMesh<GlobalOrdinal>::RendezvousMesh( const RCP_Moab& moab, 
					       const OrdinalMap& ordinal_map )
    : d_moab( moab )
    , d_ordinal_map( ordinal_map )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<typename GlobalOrdinal>
RendezvousMesh<GlobalOrdinal>::~RendezvousMesh()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Given a bounding box return the native element ordinals that are in
 * the box.
 *
 * \param box The bounding box to search with elements.
 *
 * \return The elements in the box.
 */
template<typename GlobalOrdinal>
Teuchos::Array<GlobalOrdinal> 
RendezvousMesh<GlobalOrdinal>::elementsInBox( const BoundingBox& box ) const
{
    // Get the dimension of the mesh.
    DTK_REMEMBER( moab::ErrorCode error );
    int dim = 0;
#if HAVE_DTK_DBC
    error = d_moab->get_dimension( dim );
#else
    d_moab->get_dimension( dim );
#endif
    DTK_CHECK( moab::MB_SUCCESS == error );

    // Get the mesh elements.
    moab::Range elements;
#if HAVE_DTK_DBC
    error = d_moab->get_entities_by_dimension( 0, dim, elements );
#else
    d_moab->get_entities_by_dimension( 0, dim, elements );
#endif
    DTK_CHECK( moab::MB_SUCCESS == error );
    
    // Get the elements that are in the box.
    Teuchos::Array<GlobalOrdinal> elements_in_box;
    moab::Range::iterator element_iterator;
    for ( element_iterator = elements.begin();
	  element_iterator != elements.end();
	  ++element_iterator )
    {
	if ( TopologyTools::boxElementOverlap( 
		 box, *element_iterator, d_moab ) )
	{
	    elements_in_box.push_back( getNativeOrdinal( *element_iterator ) );
	}   
    }

    return elements_in_box;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given a geometry return the native element ordinals that are in
 * the geometry.
 *
 * \param geometry The geometry to search with elements.
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
 *
 * \return The elements in the geometry.
 */
template<typename GlobalOrdinal>
template<class Geometry>
Teuchos::Array<GlobalOrdinal> 
RendezvousMesh<GlobalOrdinal>::elementsInGeometry( 
    const Geometry& geometry, const double tolerance,
    bool all_vertices_for_inclusion ) const
{
    // Get the dimension of the mesh.
    DTK_REMEMBER( moab::ErrorCode error );
    int dim = 0;
#if HAVE_DTK_DBC
    error = d_moab->get_dimension( dim );
#else
    d_moab->get_dimension( dim );
#endif
    DTK_CHECK( moab::MB_SUCCESS == error );

    // Get the mesh elements.
    moab::Range elements;
#if HAVE_DTK_DBC
    error = d_moab->get_entities_by_dimension( 0, dim, elements );
#else
    d_moab->get_entities_by_dimension( 0, dim, elements );
#endif
    DTK_CHECK( moab::MB_SUCCESS == error );
    
    // Get the elements that are in the geometry.
    Teuchos::Array<GlobalOrdinal> elements_in_geometry;
    moab::Range::iterator element_iterator;
    for ( element_iterator = elements.begin();
	  element_iterator != elements.end();
	  ++element_iterator )
    {
	if ( TopologyTools::elementInGeometry( geometry, *element_iterator, 
					       d_moab, tolerance, 
					       all_vertices_for_inclusion ) )
	{
	    elements_in_geometry.push_back( 
		getNativeOrdinal( *element_iterator ) );
	}   
    }

    return elements_in_geometry;
}

//---------------------------------------------------------------------------//
// Non-member creation methods.
//---------------------------------------------------------------------------//
/*!
 * \brief Create a RendezvousMesh from a mesh manager.
 *
 * \param mesh_manager The mesh to construct the RendezvousMesh from.
 *
 * \return The RendezvousMesh that was constructed from the mesh manager.
 */
template<typename Mesh>
Teuchos::RCP< RendezvousMesh<typename MeshTraits<Mesh>::global_ordinal_type> >
createRendezvousMeshFromMesh( const MeshManager<Mesh>& mesh_manager )
{
    // Setup types and iterators as we're outside of the class definition.
    typedef MeshTraits<Mesh> MT;
    typedef typename MT::global_ordinal_type GlobalOrdinal;
    typename MT::const_vertex_iterator vertex_iterator;
    typename MT::const_element_iterator element_iterator;

    // Setup an element handle map.
    std::map<moab::EntityHandle,GlobalOrdinal> element_ordinal_map;

    // Create a moab interface.
    DTK_REMEMBER( moab::ErrorCode error );
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );
    DTK_ENSURE( !moab.is_null() );

    // Set the mesh dimension.
    int vertex_dim = mesh_manager.dim();
#if HAVE_DTK_DBC
    error = moab->set_dimension( vertex_dim );
#else
    moab->set_dimension( vertex_dim );
#endif
    DTK_CHECK( moab::MB_SUCCESS == error );

    // Build each mesh block.
    typename MeshManager<Mesh>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	// Check the vertices and coordinates for consistency.
	GlobalOrdinal num_vertices = 
	    MeshTools<Mesh>::numVertices( *(*block_iterator) );
	DTK_REMEMBER( GlobalOrdinal num_coords = 
		       std::distance( MT::coordsBegin( *(*block_iterator) ),
				      MT::coordsEnd( *(*block_iterator) ) ) );
	DTK_CHECK( num_coords == 
		       Teuchos::as<GlobalOrdinal>(vertex_dim) * num_vertices );

	// Add the mesh vertices to moab and map the native vertex handles to
	// the moab vertex handles. This should be in a hash table. We'll need
	// one that hashes moab handles.
	double vertex_coords[3];
	Teuchos::ArrayRCP<const double> mesh_coords = 
	    MeshTools<Mesh>::coordsView( *(*block_iterator) );
	std::map<GlobalOrdinal,moab::EntityHandle> vertex_ordinal_map;
	GlobalOrdinal n = 0;
	for ( vertex_iterator = MT::verticesBegin( *(*block_iterator) );
	      vertex_iterator != MT::verticesEnd( *(*block_iterator) );
	      ++vertex_iterator, ++n )
	{
	    moab::EntityHandle moab_vertex;
	    for ( int d = 0; d < vertex_dim; ++d )
	    {
		vertex_coords[d] = mesh_coords[d*num_vertices + n];
	    }
	    for ( int d = vertex_dim; d < 3; ++d )
	    {
		vertex_coords[d] = 0.0;
	    }
#if HAVE_DTK_DBC
	    error = moab->create_vertex( vertex_coords, moab_vertex );
#else
	    moab->create_vertex( vertex_coords, moab_vertex );
#endif
	    DTK_CHECK( moab::MB_SUCCESS == error );

	    vertex_ordinal_map[ *vertex_iterator ] = moab_vertex;
	}

	// Check the elements and connectivity for consistency.
	int vertices_per_element = 
	    MT::verticesPerElement( *(*block_iterator) );
	GlobalOrdinal num_elements = 
	    MeshTools<Mesh>::numElements( *(*block_iterator) );
	DTK_REMEMBER( GlobalOrdinal num_connect = std::distance( 
			   MT::connectivityBegin( *(*block_iterator) ),
			   MT::connectivityEnd( *(*block_iterator) ) ) );
	DTK_CHECK( num_elements == num_connect / vertices_per_element &&
		       num_connect % vertices_per_element == 0 );

	// Extract the mesh elements and add them to moab.
	Teuchos::ArrayRCP<const GlobalOrdinal> mesh_connectivity = 
	    MeshTools<Mesh>::connectivityView( *(*block_iterator) );
	Teuchos::ArrayRCP<const int> permutation_list =
	    MeshTools<Mesh>::permutationView( *(*block_iterator) );
	GlobalOrdinal conn_index;
	Teuchos::Array<moab::EntityHandle> 
	    element_connectivity( vertices_per_element );

	int canonical_idx;
	n = 0;
	int element_topology = MT::elementTopology( *(*block_iterator) );
	for ( element_iterator = MT::elementsBegin( *(*block_iterator) );
	      element_iterator != MT::elementsEnd( *(*block_iterator) );
	      ++element_iterator, ++n )
	{
	    // Extract the connecting vertices for this element and apply the
	    // permutation list.
	    for ( int i = 0; i < vertices_per_element; ++i )
	    {
		canonical_idx = permutation_list[i];
		conn_index = i*num_elements + n;
		element_connectivity[ canonical_idx ] =
		    vertex_ordinal_map.find( 
			mesh_connectivity[ conn_index ] )->second;
	    }
	    DTK_CHECK( element_connectivity.size()
			   == Teuchos::as<
			       Teuchos::Array<moab::EntityHandle>::size_type>(
				   vertices_per_element) );

	    // Create the element in moab.
	    moab::EntityType entity_type = 
		moab_topology_table[ element_topology ];
	    moab::EntityHandle moab_element;
#if HAVE_DTK_DBC
	    error = moab->create_element( entity_type,
					  &element_connectivity[0],
					  element_connectivity.size(),
					  moab_element );
#else
	    moab->create_element( entity_type,
				  &element_connectivity[0],
				  element_connectivity.size(),
				  moab_element );
#endif
	    DTK_CHECK( moab::MB_SUCCESS == error );

	    // Map the moab element handle to the native element handle.
	    element_ordinal_map[ moab_element ] = *element_iterator;
	}
    }

    // Create and return the mesh.
    return Teuchos::rcp( 
	new RendezvousMesh<GlobalOrdinal>( moab, element_ordinal_map ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create a RendezvousMesh from a geometry manager. This will construct
 * a bounding region mesh from the geometric bounding boxes. Note that these
 * may overlap.
 *
 * \param geometry_manager The geometry to construct the RendezvousMesh from.
 *
 * \return The RendezvousMesh that was constructed from the geometry manager.
 */
template<typename GlobalOrdinal, typename Geometry>
Teuchos::RCP< RendezvousMesh<GlobalOrdinal> > createRendezvousMeshFromGeometry(
    const GeometryManager<Geometry,GlobalOrdinal>& geometry_manager )
{
    // Setup types and iterators as we're outside of the class definition.
    typedef GeometryTraits<Geometry> GT;

    // Setup an element handle map.
    std::map<moab::EntityHandle,GlobalOrdinal> element_ordinal_map;

    // Create a moab interface.
    DTK_REMEMBER( moab::ErrorCode error );
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );
    DTK_ENSURE( !moab.is_null() );

    // Set the mesh dimension.
    int vertex_dim = geometry_manager.dim();
#if HAVE_DTK_DBC
    error = moab->set_dimension( vertex_dim );
#else
    moab->set_dimension( vertex_dim );
#endif
    DTK_CHECK( moab::MB_SUCCESS == error );

    // Set the mesh topology.
    moab::EntityType entity_type;
    int vertices_per_element;
    if ( vertex_dim == 0 )
    {
	entity_type = moab::MBVERTEX;
	vertices_per_element = 1;
    }
    else if ( vertex_dim == 1 ) 
    {
	entity_type = moab::MBEDGE;
	vertices_per_element = 2;
    }
    else if ( vertex_dim == 2 )
    {
	entity_type = moab::MBQUAD;
	vertices_per_element = 4;
    }
    else
    {
	entity_type = moab::MBHEX;
	vertices_per_element = 8;
    }
    
    // Extract the geometry bounding boxes from the manager. We will turn
    // these into hexahedrons, squares, or lines.
    std::map<GlobalOrdinal,moab::EntityHandle> vertex_ordinal_map;
    Teuchos::Tuple<double,6> geom_bounds;
    Teuchos::Array<double> vertex_coords(3*vertices_per_element);
    Teuchos::ArrayRCP<Geometry> geometry = geometry_manager.geometry();
    typename Teuchos::ArrayRCP<Geometry>::const_iterator geometry_begin = 
	geometry.begin();
    typename Teuchos::ArrayRCP<Geometry>::const_iterator geometry_iterator;
    for ( geometry_iterator = geometry.begin(); 
	  geometry_iterator != geometry.end();
	  ++geometry_iterator )
    {
	// Get the bounding coordinates.
	geom_bounds = GT::boundingBox( *geometry_iterator ).getBounds();
	
	// 0D vertex case.
	if ( vertex_dim == 0 )
	{
	    vertex_coords[0] = geom_bounds[0];
	}
	// 1D Line case.
	else if ( vertex_dim == 1 )
	{
	    vertex_coords[0] = geom_bounds[0];
	    vertex_coords[1] = 0.0;
	    vertex_coords[2] = 0.0;

	    vertex_coords[3] = geom_bounds[3];
	    vertex_coords[4] = 0.0;
	    vertex_coords[5] = 0.0;
	}
	// 2D Quad case.
	else if ( vertex_dim == 2 )
	{
	    vertex_coords[0] = geom_bounds[0];
	    vertex_coords[1] = geom_bounds[1];
	    vertex_coords[2] = 0.0;

	    vertex_coords[3] = geom_bounds[3];
	    vertex_coords[4] = geom_bounds[1];
	    vertex_coords[5] = 0.0;

	    vertex_coords[6] = geom_bounds[3];
	    vertex_coords[7] = geom_bounds[4];
	    vertex_coords[8] = 0.0;

	    vertex_coords[9]  = geom_bounds[0];
	    vertex_coords[10] = geom_bounds[4];
	    vertex_coords[11] = 0.0;
	}
	// 3D hex case.
	else 
	{
	    vertex_coords[0] = geom_bounds[0];
	    vertex_coords[1] = geom_bounds[1];
	    vertex_coords[2] = geom_bounds[2];

	    vertex_coords[3] = geom_bounds[3];
	    vertex_coords[4] = geom_bounds[1];
	    vertex_coords[5] = geom_bounds[2];

	    vertex_coords[6] = geom_bounds[3];
	    vertex_coords[7] = geom_bounds[4];
	    vertex_coords[8] = geom_bounds[2];

	    vertex_coords[9]  = geom_bounds[0];
	    vertex_coords[10] = geom_bounds[4];
	    vertex_coords[11] = geom_bounds[2];

	    vertex_coords[12] = geom_bounds[0];
	    vertex_coords[13] = geom_bounds[1];
	    vertex_coords[14] = geom_bounds[5];

	    vertex_coords[15] = geom_bounds[3];
	    vertex_coords[16] = geom_bounds[1];
	    vertex_coords[17] = geom_bounds[5];

	    vertex_coords[18] = geom_bounds[3];
	    vertex_coords[19] = geom_bounds[4];
	    vertex_coords[20] = geom_bounds[5];

	    vertex_coords[21] = geom_bounds[0];
	    vertex_coords[22] = geom_bounds[4];
	    vertex_coords[23] = geom_bounds[5];
	}

	// Build the vertices.
	moab::Range moab_vertices;
#if HAVE_DTK_DBC
	error = moab->create_vertices( &vertex_coords[0], vertices_per_element,
				       moab_vertices );
#else
	moab->create_vertices( &vertex_coords[0], vertices_per_element,
			       moab_vertices );
#endif
	DTK_CHECK( moab::MB_SUCCESS == error );
	DTK_CHECK( (int) moab_vertices.size() == vertices_per_element );
	Teuchos::Array<moab::EntityHandle> vert_copy( vertices_per_element );
	std::copy( moab_vertices.begin(), moab_vertices.end(), 
		   vert_copy.begin() );
	moab_vertices.clear();
	
	// Build the bounding entity.
	moab::EntityHandle moab_element;
#if HAVE_DTK_DBC
	error = moab->create_element( entity_type,
				      &vert_copy[0],
				      vertices_per_element,
				      moab_element );
#else
	moab->create_element( entity_type,
			      &vert_copy[0],
			      vertices_per_element,
			      moab_element );
#endif
	DTK_CHECK( moab::MB_SUCCESS == error );

	// Map the moab element handle to the local geometry ordinal.
	element_ordinal_map[ moab_element ] = 
	    std::distance( geometry_begin, geometry_iterator );
    }

    // Create and return the mesh.
    return Teuchos::rcp( 
	new RendezvousMesh<GlobalOrdinal>( moab, element_ordinal_map ) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_RENDEZVOUSMESH_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_RendezvousMesh_def.hpp
//---------------------------------------------------------------------------//

