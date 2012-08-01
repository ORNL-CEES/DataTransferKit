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
#include <DTK_Exception.hpp>

#include <MBCore.hpp>

#include <Teuchos_ENull.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
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
 */
template<class Mesh>
Teuchos::RCP< RendezvousMesh<typename MeshTraits<Mesh>::global_ordinal_type> >
createRendezvousMesh( const MeshManager<Mesh>& mesh_manager )
{
    // Setup types and iterators as we're outside of the class definition.
    typedef MeshTraits<Mesh> MT;
    typedef typename MT::global_ordinal_type GlobalOrdinal;
    typename MT::const_node_iterator node_iterator;
    typename MT::const_element_iterator element_iterator;

    // Setup an element handle map.
    std::map<moab::EntityHandle,GlobalOrdinal> element_handle_map;

    // Create a moab interface.
    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );
    testPostcondition( moab != Teuchos::null,
		       "Error creating MOAB interface" );

    // Set the mesh dimension.
    std::size_t node_dim = mesh_manager.dim();
    error = moab->set_dimension( node_dim );
    testInvariant( moab::MB_SUCCESS == error, 
		   "Failed to set MOAB mesh dimension: " +
		   moab->get_error_string( error ) );

    // Build each mesh block.
    typename MeshManager<Mesh>::BlockIterator block_iterator;
    for ( block_iterator = mesh_manager.blocksBegin();
	  block_iterator != mesh_manager.blocksEnd();
	  ++block_iterator )
    {
	// Check the nodes and coordinates for consistency.
	GlobalOrdinal num_nodes = MeshTools<Mesh>::numNodes( *block_iterator );
	GlobalOrdinal num_coords = 
	    std::distance( MT::coordsBegin( *block_iterator ),
			   MT::coordsEnd( *block_iterator ) );
	testInvariant( 
	    num_coords == (GlobalOrdinal) node_dim * num_nodes,
	    "Number of coordinates provided != node_dim * number of nodes" );

	// Add the mesh nodes to moab and map the native vertex handles to the
	// moab vertex handles. This should be in a hash table. We'll need one
	// that hashes moab handles.
	double vertex_coords[3];
	Teuchos::ArrayRCP<const double> mesh_coords = 
	    MeshTools<Mesh>::coordsView( *block_iterator );
	std::map<GlobalOrdinal,moab::EntityHandle> vertex_handle_map;
	GlobalOrdinal n = 0;
	for ( node_iterator = MT::nodesBegin( *block_iterator );
	      node_iterator != MT::nodesEnd( *block_iterator );
	      ++node_iterator, ++n )
	{
	    moab::EntityHandle moab_vertex;
	    for ( std::size_t d = 0; d < node_dim; ++d )
	    {
		vertex_coords[d] = mesh_coords[d*num_nodes + n];
	    }
	    for ( std::size_t d = node_dim; d < 3; ++d )
	    {
		vertex_coords[d] = 0.0;
	    }
	    error = moab->create_vertex( vertex_coords, moab_vertex );
	    testInvariant( moab::MB_SUCCESS == error, 
			   "Failed to create vertices in MOAB: " +
			   moab->get_error_string( error ) );
	    vertex_handle_map[ *node_iterator ] = moab_vertex;
	}

	// Check the elements and connectivity for consistency.
	int nodes_per_element = 
	    MT::nodesPerElement( *block_iterator );
	GlobalOrdinal num_elements = 
	    MeshTools<Mesh>::numElements( *block_iterator );
	GlobalOrdinal num_connect = 
	    std::distance( MT::connectivityBegin( *block_iterator ),
			   MT::connectivityEnd( *block_iterator ) );
	testPrecondition( 
	    num_elements == num_connect / nodes_per_element &&
	    num_connect % nodes_per_element == 0,
	    "Connectivity array inconsistent with element description." );

	// Extract the mesh elements and add them to moab.
	Teuchos::ArrayRCP<const GlobalOrdinal> mesh_connectivity = 
	    MeshTools<Mesh>::connectivityView( *block_iterator );
	Teuchos::ArrayRCP<const std::size_t> permutation_list =
	    MeshTools<Mesh>::permutationView( *block_iterator );
	GlobalOrdinal conn_index;
	Teuchos::Array<moab::EntityHandle> 
	    element_connectivity( nodes_per_element );

	std::size_t canonical_idx;
	n = 0;
	int element_topology = MT::elementTopology( *block_iterator );
	for ( element_iterator = MT::elementsBegin( *block_iterator );
	      element_iterator != MT::elementsEnd( *block_iterator );
	      ++element_iterator, ++n )
	{
	    // Extract the connecting nodes for this element and apply the
	    // permutation list.
	    for ( int i = 0; i < nodes_per_element; ++i )
	    {
		canonical_idx = permutation_list[i];
		conn_index = i*num_elements + n;
		element_connectivity[ canonical_idx ] =
		    vertex_handle_map.find( 
			mesh_connectivity[ conn_index ] )->second;
	    }
	    testInvariant( (int) element_connectivity.size() == nodes_per_element,
			   "Element connectivity size != nodes per element." );

	    // Create the element in moab.
	    moab::EntityType entity_type = 
		moab_topology_table[ element_topology ];
	    moab::EntityHandle moab_element;
	    error = moab->create_element( entity_type,
					  &element_connectivity[0],
					  element_connectivity.size(),
					  moab_element );
	    testInvariant( moab::MB_SUCCESS == error,
			   "Failed to create element in MOAB:" +
			   moab->get_error_string( error ) );

	    // Map the moab element handle to the native element handle.
	    element_handle_map[ moab_element ] = *element_iterator;
	}
    }

    // Temporary debug output.
    moab::Range elements;
    error = moab->get_entities_by_dimension( 0, 3, elements );
    testInvariant( moab::MB_SUCCESS == error,
		   "Failed to get mesh elements: " + 
		   moab->get_error_string( error ) );

    error = moab->list_entities( elements );
    testInvariant( moab::MB_SUCCESS == error,
		   "Failed to get mesh elements: " + 
		   moab->get_error_string( error ) );

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

