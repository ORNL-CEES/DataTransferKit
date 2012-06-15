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
#include <DTK_Exception.hpp>

#include <MBCore.hpp>

#include <Teuchos_ENull.hpp>
#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename GlobalOrdinal>
RendezvousMesh<GlobalOrdinal>::RendezvousMesh( const RCP_Moab& moab, 
					       const moab::Range& elements,
					       const HandleMap& handle_map )
    : d_moab( moab )
    , d_elements( elements )
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
 * \brief Create a RendezvousMesh from an object that implements mesh traits.
 */
template<class Mesh>
Teuchos::RCP< RendezvousMesh<typename MeshTraits<Mesh>::global_ordinal_type> >
createRendezvousMesh( const Mesh& mesh )
{
    // Setup types and iterators as we're outside of the class definition.
    typedef MeshTraits<Mesh> MT;
    typedef typename MT::global_ordinal_type GlobalOrdinal;
    typename MT::const_node_iterator node_iterator;
    typename MT::const_coordinate_iterator coord_iterator;
    typename MT::const_element_iterator element_iterator;
    typename MT::const_connectivity_iterator conn_iterator;

    // Create a moab interface.
    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );
    testPostcondition( moab != Teuchos::null,
		       "Error creating MOAB interface" );

    // Check the nodes and coordinates for consistency.
    std::size_t node_dim = MT::nodeDim( mesh );
    GlobalOrdinal num_nodes = std::distance( MT::nodesBegin( mesh ), 
					     MT::nodesEnd( mesh ) );
    GlobalOrdinal num_coords = std::distance( MT::coordsBegin( mesh ),
					      MT::coordsEnd( mesh ) );
    testInvariant( 
	num_coords == (GlobalOrdinal) node_dim * num_nodes,
	"Number of coordinates provided != node_dim * number of nodes" );

    // Add the mesh nodes to moab and map the native vertex handles to the
    // moab vertex handles. This should be in a hash table. We'll need one
    // that hashes moab handles.
    double vertex_coords[3];
    coord_iterator = MT::coordsBegin( mesh );
    std::map<GlobalOrdinal,moab::EntityHandle> vertex_handle_map;
    GlobalOrdinal n = 0;
    for ( node_iterator = MT::nodesBegin( mesh );
	  node_iterator != MT::nodesEnd( mesh );
	  ++node_iterator, ++n )
    {
	moab::EntityHandle moab_vertex;
	for ( std::size_t d = 0; d < node_dim; ++d )
	{
	    vertex_coords[d] = coord_iterator[d*num_nodes + n];
	}
	for ( std::size_t d = node_dim; d < 3; ++d )
	{
	    vertex_coords[d] = 0.0;
	}
	error = moab->create_vertex( vertex_coords, moab_vertex );
	testInvariant( moab::MB_SUCCESS == error, 
		       "Failed to create vertices in MOAB." );
	vertex_handle_map[ *node_iterator ] = moab_vertex;
    }

    // Check the elements and connectivity for consistency.
    int nodes_per_element = 
	MT::nodesPerElement( mesh );
    GlobalOrdinal num_elements = std::distance( MT::elementsBegin( mesh ),
						MT::elementsEnd( mesh ) );
    GlobalOrdinal num_connect = std::distance( MT::connectivityBegin( mesh ),
					       MT::connectivityEnd( mesh ) );
    testPrecondition( 
	num_elements == num_connect / nodes_per_element &&
	num_connect % nodes_per_element == 0,
	"Connectivity array inconsistent with element description." );

    // Extract the mesh elements and add them to moab.
    conn_iterator = MT::connectivityBegin( mesh );
    GlobalOrdinal conn_index;
    moab::Range moab_elements;
    Teuchos::Array<moab::EntityHandle> element_connectivity;
    std::map<moab::EntityHandle,GlobalOrdinal> element_handle_map;
    n = 0;
    for ( element_iterator = MT::elementsBegin( mesh );
	  element_iterator != MT::elementsEnd( mesh );
	  ++element_iterator, ++n )
    {
	// Extract the connecting nodes for this element.
	element_connectivity.clear();
	for ( int i = 0; i < nodes_per_element; ++i )
	{
	    conn_index = i*num_elements + n;
	    element_connectivity.push_back( 
		vertex_handle_map.find( conn_iterator[ conn_index ] )->second );
	}
	testInvariant( (int) element_connectivity.size() == nodes_per_element,
		       "Element connectivity size != nodes per element." );

	// Create the element in moab.
	moab::EntityType entity_type = moab_topology_table[ 
	    MT::elementTopology( mesh ) ];
	moab::EntityHandle moab_element;
	error = moab->create_element( entity_type,
				      &element_connectivity[0],
				      element_connectivity.size(),
				      moab_element );
	testInvariant( moab::MB_SUCCESS == error,
		       "Failed to create element in MOAB." );
	moab_elements.insert( moab_element );

	// Map the moab element handle to the native element handle.
	element_handle_map[ moab_element ] = *element_iterator;
    }
    
    // Create and return the mesh.
    return Teuchos::rcp( new RendezvousMesh<GlobalOrdinal>( 
			     moab, moab_elements, element_handle_map ) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_RENDEZVOUSMESH_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_RendezvousMesh_def.hpp
//---------------------------------------------------------------------------//

