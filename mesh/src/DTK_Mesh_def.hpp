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
#include <DTK_NodeTraits.hpp>
#include <DTK_ElementTraits.hpp>
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
// Non-member creation functions.
//---------------------------------------------------------------------------//
/*!
 * \brief Create a mesh from a DataSource.
 */
template<typename NodeField, typename ElementField, typename DataField>
Teuchos::RCP< Mesh<typename ElementField::value_type::handle_type> >
createMeshFromDataSource(
    const Teuchos::RCP< 
	DataSource<NodeField,ElementField,DataField> >& data_source )
{
    // Setup for data source.
    typedef DataSource<NodeField,ElementField,DataField> DS;

    // Setup for node types.
    typedef typename DS::node_type node_type;
    typedef typename DS::node_handle_type node_handle_type;
    typedef typename DS::node_coordinate_type coordinate_type;
    typename FieldTraits<NodeField>::const_iterator node_iterator;    
    typename NodeTraits<node_type>::const_coordinate_iterator 
	coord_iterator;

    // Setup for element types.
    typedef typename DS::element_type element_type;
    typedef typename DS::element_handle_type element_handle_type;
    typename FieldTraits<ElementField>::const_iterator element_iterator;
    typename ElementTraits<element_type>::const_connectivity_iterator
	conn_iterator;

    // Create a moab interface.
    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );
    testPostcondition( moab != Teuchos::null,
		       "Error creating MOAB interface" );

    // Extract the source mesh nodes;
    int zeros_to_add = 3 - NodeTraits<node_type>::dim();
    std::vector<node_handle_type> node_handles;
    std::vector<coordinate_type> node_coords;
    NodeField source_nodes = data_source->getSourceMeshNodes();
    for ( node_iterator = FieldTraits<NodeField>::begin( source_nodes );
	  node_iterator != FieldTraits<NodeField>::end( source_nodes );
	  ++node_iterator )
    {
	node_handles.push_back( 
	    NodeTraits<node_type>::handle( *node_iterator ) );

	for ( coord_iterator = 
		  NodeTraits<node_type>::coordsBegin( *node_iterator );
	      coord_iterator != 
		  NodeTraits<node_type>::coordsEnd( *node_iterator );
	      ++coord_iterator )
	{
	    node_coords.push_back( *coord_iterator );
	}

	for ( int i = 0; i < zeros_to_add; ++i )
	{
	    node_coords.push_back( 0.0 );
	}
    }
    assert( node_coords.size() == 
	    3 * FieldTraits<NodeField>::size( source_nodes ) );

    // Add the source mesh nodes to moab.
    moab::Range vertices;
    error = moab->create_vertices( &node_coords[0], node_handles.size(),
				   vertices );
    testInvariant( moab::MB_SUCCESS == error, 
		   "Failed to create vertices in MOAB." );
    testPostcondition( !vertices.empty(),
		       "Vertex range is empty." );
    assert( vertices.size() == node_handles.size() );

    // Map the native vertex handles to the moab vertex handles.
    moab::Range::const_iterator range_iterator;
    typename std::vector<node_handle_type>::const_iterator handle_iterator;
    std::map<node_handle_type,moab::EntityHandle> vertex_handle_map;
    for ( range_iterator = vertices.begin(),
	 handle_iterator = node_handles.begin();
	  range_iterator != vertices.end();
	  ++range_iterator, ++handle_iterator )
    {
	vertex_handle_map[ *handle_iterator ] = *range_iterator;
    }

    // Extract the source mesh elements and add them to moab.
    moab::Range elements;
    std::vector<moab::EntityHandle> element_connectivity;
    std::map<moab::EntityHandle,element_handle_type> element_handle_map;
    ElementField source_elements = data_source->getSourceMeshElements();
    for ( element_iterator = FieldTraits<ElementField>::begin( source_elements );
	  element_iterator != FieldTraits<ElementField>::end( source_elements );
	  ++element_iterator )
    {
	element_connectivity.clear();
	for ( conn_iterator = 
		  ElementTraits<element_type>::connectivityBegin( *element_iterator );
	      conn_iterator != 
		  ElementTraits<element_type>::connectivityEnd( *element_iterator );
	      ++conn_iterator )
	{
	    element_connectivity.push_back( vertex_handle_map[*conn_iterator] );
	}

	testInvariant( element_connectivity.size() == 
		       ElementTraits<element_type>::numNodes(),
		       "Element connectivity size != number of element nodes." );

	moab::EntityType entity_type = moab_topology_table[ 
	    ElementTraits<element_type>::topology() ];
	moab::EntityHandle moab_element;
	error = moab->create_element( entity_type,
				      &element_connectivity[0],
				      element_connectivity.size(),
				      moab_element );
	testInvariant( moab::MB_SUCCESS == error,
		       "Failed to create element in MOAB." );

	elements.insert( moab_element );

	element_handle_map[ moab_element ] =
	    ElementTraits<element_type>::handle( *element_iterator );
    }
    
    // Create and return the mesh.
    return Teuchos::rcp( 
	new Mesh<element_handle_type>( moab, elements, element_handle_map ) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MESH_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_Mesh_def.hpp
//---------------------------------------------------------------------------//

