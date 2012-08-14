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
#include <cassert>

#include "DTK_MeshTraits.hpp"
#include "DTK_MeshTools.hpp"
#include "DTK_Assertion.hpp"

#include <MBCore.hpp>

#include <Teuchos_ENull.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 * 
 * \param moab The Moab interface to build the RendezvousMesh with.
 *
 * \param handle_map A map relating the client element global ordinals to the
 * Moab element handles.
 */
template<typename GlobalOrdinal>
RendezvousMesh<GlobalOrdinal>::RendezvousMesh( const RCP_Moab& moab, 
					       const HandleMap& handle_map )
    : d_moab( moab )
    , d_handle_map( handle_map )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<typename GlobalOrdinal>
RendezvousMesh<GlobalOrdinal>::~RendezvousMesh()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Non-member creation methods.
//---------------------------------------------------------------------------//
/*!
 * \brief Create a RendezvousMesh from a mesh manager.
 *
 * \param mesh_manager The mesh to construct the RendezvousMesh from.
 *
 * \return The RendezvousMesh that was constructed from the mesh.
 */
template<class Mesh>
Teuchos::RCP< RendezvousMesh<typename MeshTraits<Mesh>::global_ordinal_type> >
createRendezvousMesh( const MeshManager<Mesh>& mesh_manager )
{
    // Setup types and iterators as we're outside of the class definition.
    typedef MeshTraits<Mesh> MT;
    typedef typename MT::global_ordinal_type GlobalOrdinal;
    typename MT::const_vertex_iterator vertex_iterator;
    typename MT::const_element_iterator element_iterator;

    // Setup an element handle map.
    std::map<moab::EntityHandle,GlobalOrdinal> element_handle_map;

    // Create a moab interface.
    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );
    testPostcondition( moab != Teuchos::null );

    // Set the mesh dimension.
    int vertex_dim = mesh_manager.dim();
    error = moab->set_dimension( vertex_dim );
    assert( moab::MB_SUCCESS == error );

    // Build each mesh block.
    typename MeshManager<Mesh>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	// Check the vertices and coordinates for consistency.
	GlobalOrdinal num_vertices = 
	    MeshTools<Mesh>::numVertices( *block_iterator );
	rememberValue( GlobalOrdinal num_coords = 
		       std::distance( MT::coordsBegin( *block_iterator ),
				      MT::coordsEnd( *block_iterator ) ) );
	testInvariant( num_coords == 
		       (GlobalOrdinal) vertex_dim * num_vertices );

	// Add the mesh vertices to moab and map the native vertex handles to
	// the moab vertex handles. This should be in a hash table. We'll need
	// one that hashes moab handles.
	double vertex_coords[3];
	Teuchos::ArrayRCP<const double> mesh_coords = 
	    MeshTools<Mesh>::coordsView( *block_iterator );
	std::map<GlobalOrdinal,moab::EntityHandle> vertex_handle_map;
	GlobalOrdinal n = 0;
	for ( vertex_iterator = MT::verticesBegin( *block_iterator );
	      vertex_iterator != MT::verticesEnd( *block_iterator );
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
	    error = moab->create_vertex( vertex_coords, moab_vertex );
	    assert( moab::MB_SUCCESS == error );

	    vertex_handle_map[ *vertex_iterator ] = moab_vertex;
	}

	// Check the elements and connectivity for consistency.
	int vertices_per_element = 
	    MT::verticesPerElement( *block_iterator );
	GlobalOrdinal num_elements = 
	    MeshTools<Mesh>::numElements( *block_iterator );
	rememberValue( GlobalOrdinal num_connect = 
		       std::distance( MT::connectivityBegin( *block_iterator ),
				      MT::connectivityEnd( *block_iterator ) ) );
	testInvariant( num_elements == num_connect / vertices_per_element &&
		       num_connect % vertices_per_element == 0 );

	// Extract the mesh elements and add them to moab.
	Teuchos::ArrayRCP<const GlobalOrdinal> mesh_connectivity = 
	    MeshTools<Mesh>::connectivityView( *block_iterator );
	Teuchos::ArrayRCP<const int> permutation_list =
	    MeshTools<Mesh>::permutationView( *block_iterator );
	GlobalOrdinal conn_index;
	Teuchos::Array<moab::EntityHandle> 
	    element_connectivity( vertices_per_element );

	int canonical_idx;
	n = 0;
	int element_topology = MT::elementTopology( *block_iterator );
	for ( element_iterator = MT::elementsBegin( *block_iterator );
	      element_iterator != MT::elementsEnd( *block_iterator );
	      ++element_iterator, ++n )
	{
	    // Extract the connecting vertices for this element and apply the
	    // permutation list.
	    for ( int i = 0; i < vertices_per_element; ++i )
	    {
		canonical_idx = permutation_list[i];
		conn_index = i*num_elements + n;
		element_connectivity[ canonical_idx ] =
		    vertex_handle_map.find( 
			mesh_connectivity[ conn_index ] )->second;
	    }
	    testInvariant( (int) element_connectivity.size() 
			   == vertices_per_element );

	    // Create the element in moab.
	    moab::EntityType entity_type = 
		moab_topology_table[ element_topology ];
	    moab::EntityHandle moab_element;
	    error = moab->create_element( entity_type,
					  &element_connectivity[0],
					  element_connectivity.size(),
					  moab_element );
	    assert( moab::MB_SUCCESS == error );

	    // Map the moab element handle to the native element handle.
	    element_handle_map[ moab_element ] = *element_iterator;
	}
    }

    // Create and return the mesh.
    return Teuchos::rcp( 
	new RendezvousMesh<GlobalOrdinal>( moab, element_handle_map ) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_RENDEZVOUSMESH_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_RendezvousMesh_def.hpp
//---------------------------------------------------------------------------//

