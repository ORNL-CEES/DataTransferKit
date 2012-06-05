//---------------------------------------------------------------------------//
/*!
 * \file DTK_Mesh_def.hpp
 * \author Stuart R. Slattery
 * \brief Concrete mesh template definitions.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESH_DEF_HPP
#define DTK_MESH_DEF_HPP

#include <vector>
#include <cassert>

#include <DTK_Exception.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_FieldTraits.hpp>

#include <MBCore.hpp>

#include <Teuchos_ENull.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename ElementHandle>
Mesh<ElementHandle>::Mesh( const RCP_Moab& moab, 
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
template<typename ElementHandle>
Mesh<ElementHandle>::~Mesh()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Non-member creation methods.
//---------------------------------------------------------------------------//
/*!
 * \brief Create a mesh from a node and element field.
 */
template<typename MeshObject>
Teuchos::RCP< Mesh<typename MeshObject::handle_type> > 
createMesh( const MeshObject& mesh_object )
{
    // Setup types and iterators
    typedef typename MeshTraits<MeshObject>::handle_type handle_type;
    typename MeshTraits<MeshObject>::const_handle_iterator handle_iterator;
    typename MeshTraits<MeshObject>::const_handle_iterator conn_iterator;
    typename MeshTraits<MeshObject>::const_coordinate_iterator coord_iterator;

    // Create a moab interface.
    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );
    testPostcondition( moab != Teuchos::null,
		       "Error creating MOAB interface" );


    // Check the nodes and coordinates for consistency.
    int num_nodes = 
	std::distance( MeshTraits<MeshObject>::nodesBegin( mesh_object ),
		       MeshTraits<MeshObject>::nodesEnd( mesh_object ) );
    int num_coords = 
	std::distance( MeshTraits<MeshObject>::coordsBegin( mesh_object ),
		       MeshTraits<MeshObject>::coordsEnd( mesh_object ) );
    testInvariant( num_coords == 3 * num_nodes,
		   "Number of coordinates provided != 3 * number of nodes" );

    // Add the source mesh nodes to moab.    
    moab::Range vertices;
    error = moab->create_vertices(
	&( *MeshTraits<MeshObject>::coordsBegin( mesh_object ) ),
	num_nodes, vertices );
    testInvariant( moab::MB_SUCCESS == error, 
		   "Failed to create vertices in MOAB." );
    testPostcondition( !vertices.empty(),
		       "Vertex range is empty." );
    assert( (int) vertices.size() == num_nodes );

    // Map the native vertex handles to the moab vertex handles.
    moab::Range::const_iterator range_iterator;
    std::map<handle_type,moab::EntityHandle> vertex_handle_map;
    for ( range_iterator = vertices.begin(),
	 handle_iterator = MeshTraits<MeshObject>::nodesBegin( mesh_object );
	  range_iterator != vertices.end();
	  ++range_iterator, ++handle_iterator )
    {
	vertex_handle_map[ *handle_iterator ] = *range_iterator;
    }

    // Check the elements and connectivity for consistency.
    int nodes_per_element = 
	MeshTraits<MeshObject>::nodesPerElement( mesh_object );
    int num_elements = 
	std::distance( MeshTraits<MeshObject>::elementsBegin( mesh_object ),
		       MeshTraits<MeshObject>::elementsEnd( mesh_object ) );
    int num_connect =
	std::distance( MeshTraits<MeshObject>::connectivityBegin( mesh_object ),
		       MeshTraits<MeshObject>::connectivityEnd( mesh_object ) );
    testInvariant( num_elements == num_connect / nodes_per_element &&
		   num_connect % nodes_per_element == 0,
		   "Connectivity array inconsistent with element description." );

    // Extract the source mesh elements and add them to moab.
    moab::Range moab_elements;
    std::vector<moab::EntityHandle> element_connectivity;
    std::map<moab::EntityHandle,handle_type> element_handle_map;
    for ( handle_iterator = MeshTraits<MeshObject>::elementsBegin( mesh_object ),
	    conn_iterator = MeshTraits<MeshObject>::connectivityBegin( mesh_object );
	  handle_iterator != MeshTraits<MeshObject>::elementsEnd( mesh_object );
	  ++handle_iterator )
    {
	// Extract the connecting nodes for this element.
	element_connectivity.clear();
	for ( int n = 0; n < nodes_per_element; ++n, ++conn_iterator )
	{
	    element_connectivity.push_back( vertex_handle_map[*conn_iterator] );
	}
	testInvariant( (int) element_connectivity.size() == nodes_per_element,
		       "Element connectivity size != nodes per element." );

	// Creat the element in moab.
	moab::EntityType entity_type = moab_topology_table[ 
	    MeshTraits<MeshObject>::elementTopology( mesh_object ) ];
	moab::EntityHandle moab_element;
	error = moab->create_element( entity_type,
				      &element_connectivity[0],
				      element_connectivity.size(),
				      moab_element );
	testInvariant( moab::MB_SUCCESS == error,
		       "Failed to create element in MOAB." );
	moab_elements.insert( moab_element );

	// Map the moab element handle to the native element handle.
	element_handle_map[ moab_element ] = *handle_iterator;
    }
    
    // Create and return the mesh.
    return Teuchos::rcp( 
	new Mesh<handle_type>( moab, moab_elements, element_handle_map ) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MESH_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_Mesh_def.hpp
//---------------------------------------------------------------------------//

