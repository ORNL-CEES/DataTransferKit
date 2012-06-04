//---------------------------------------------------------------------------//
/*!
 * \file tstMesh.cpp
 * \author Stuart R. Slattery
 * \brief Mesh unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_Mesh.hpp>
#include <DTK_DataSource.hpp>
#include <DTK_CoreTypes.hpp>
#include <DTK_NodeTraits.hpp>
#include <DTK_ElementTraits.hpp>
#include <DTK_FieldTraits.hpp>
#include <DTK_Exception.hpp>

#include <mpi.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>

#include <MBRange.hpp>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Node Implementation
//---------------------------------------------------------------------------//

class MyNode
{
  private:

    std::size_t d_handle;
    std::vector<double> d_coords;

  public:

    typedef int    handle_type;
    typedef double coordinate_type;
    
    MyNode( double x, double y, double z, int handle )
	: d_handle( handle )
    {
	d_coords.push_back(x);
	d_coords.push_back(y);
	d_coords.push_back(z);
    }

    ~MyNode()
    { /* ... */ }

    int handle() const
    { return d_handle; }

    std::vector<double>::const_iterator coordsBegin() const
    { return d_coords.begin(); }

    std::vector<double>::const_iterator coordsEnd() const
    { return d_coords.end(); }
};

//---------------------------------------------------------------------------//
// Element Implementation
//---------------------------------------------------------------------------//

class MyHex
{
  private:

    std::size_t d_handle;
    std::vector<int> d_connectivity;

  public:

    typedef std::size_t handle_type;

    MyHex( int node_0, int node_1, int node_2, int node_3,
	   int node_4, int node_5, int node_6, int node_7,
	   std::size_t handle )
	: d_handle( handle )
    {
	d_connectivity.push_back( node_0 );
	d_connectivity.push_back( node_1 );
	d_connectivity.push_back( node_2 );
	d_connectivity.push_back( node_3 );
	d_connectivity.push_back( node_4 );
	d_connectivity.push_back( node_5 );
	d_connectivity.push_back( node_6 );
	d_connectivity.push_back( node_7 );
    }

    ~MyHex()
    { /* ... */ }

    int handle() const
    { return d_handle; }

    std::vector<int>::const_iterator connectivityBegin() const
    { return d_connectivity.begin(); }

    std::vector<int>::const_iterator connectivityEnd() const
    { return d_connectivity.end(); }
};

//---------------------------------------------------------------------------//
// DTK Traits Specializations
//---------------------------------------------------------------------------//
namespace DataTransferKit
{

//---------------------------------------------------------------------------//
// NodeTraits specialization for the MyNode implementation.
template<>
struct NodeTraits<MyNode>
{
    typedef typename MyNode::handle_type                 handle_type;
    typedef typename MyNode::coordinate_type             coordinate_type;
    typedef typename std::vector<double>::const_iterator 
    const_coordinate_iterator;
    
    static inline std::size_t dim()
    { return 3;}
    
    static inline handle_type handle( const MyNode& node ) 
    { return node.handle(); }
    
    static inline const_coordinate_iterator coordsBegin( const MyNode& node ) 
    { return node.coordsBegin(); }

    static inline const_coordinate_iterator coordsEnd( const MyNode& node ) 
    { return node.coordsEnd(); }
};

//---------------------------------------------------------------------------//
// ElementTraits specialization for the MyHex implementation.
template<>
struct ElementTraits<MyHex>
{
    typedef typename MyHex::handle_type              handle_type;
    typedef typename std::vector<int>::const_iterator 
    const_connectivity_iterator;

    static inline std::size_t type()
    { return DTK_REGION; }

    static inline std::size_t topology()
    { return DTK_HEXAHEDRON; }

    static inline std::size_t numNodes()
    { return 8; }

    static inline handle_type handle( const MyHex &hex )
    { return hex.handle(); }

    static inline const_connectivity_iterator 
    connectivityBegin( const MyHex &hex )
    { return hex.connectivityBegin(); }

    static inline const_connectivity_iterator 
    connectivityEnd( const MyHex &hex )
    { return hex.connectivityEnd(); }
};

//---------------------------------------------------------------------------//
// FieldTraits specialization for the node field.
template<>
struct FieldTraits< std::vector<MyNode> >
{
    typedef MyNode                                value_type;
    typedef std::vector<MyNode>::iterator         iterator;
    typedef std::vector<MyNode>::const_iterator   const_iterator;
    
    static inline std::size_t size( const std::vector<MyNode> &node_field )
    { return node_field.size(); }

    static iterator begin( std::vector<MyNode> &node_field )
    { return node_field.begin(); }

    static const_iterator begin( const std::vector<MyNode> &node_field )
    { return node_field.begin(); }

    static inline iterator end( std::vector<MyNode> &node_field )
    { return node_field.end(); }

    static inline const_iterator end( const std::vector<MyNode> &node_field )
    { return node_field.end(); }

    static inline bool empty( const std::vector<MyNode> &node_field )
    { return node_field.empty(); }
};

//---------------------------------------------------------------------------//
// FieldTraits specialization for the element field.
template<>
struct FieldTraits< std::vector<MyHex> >
{
    typedef MyHex                                value_type;
    typedef std::vector<MyHex>::iterator         iterator;
    typedef std::vector<MyHex>::const_iterator   const_iterator;
    
    static inline std::size_t size( const std::vector<MyHex> &quad_field )
    { return quad_field.size(); }

    static inline iterator begin( std::vector<MyHex> &quad_field )
    { return quad_field.begin(); }

    static inline const_iterator begin( const std::vector<MyHex> &quad_field )
    { return quad_field.begin(); }

    static inline iterator end( std::vector<MyHex> &quad_field )
    { return quad_field.end(); }

    static inline const_iterator end( const std::vector<MyHex> &quad_field )
    { return quad_field.end(); }

    static inline bool empty(  const std::vector<MyHex> &quad_field )
    { return quad_field.empty(); }
};

//---------------------------------------------------------------------------//
// FieldTraits specialization for the data field.
template<>
struct FieldTraits< std::vector<double> >
{
    typedef MyNode value_type;
    typedef std::vector<double>::iterator iterator;
    typedef std::vector<double>::const_iterator const_iterator;
    
    static inline std::size_t size( const std::vector<double> &data_field )
    { return data_field.size(); }

    static inline iterator begin( std::vector<double> &data_field )
    { return data_field.begin(); }

    static inline const_iterator begin( const std::vector<double> &data_field )
    { return data_field.begin(); }

    static inline iterator end( std::vector<double> &data_field )
    { return data_field.end(); }

    static inline const_iterator end( const std::vector<double> &data_field )
    { return data_field.end(); }

    static inline bool empty( const std::vector<double> &data_field )
    { return data_field.empty(); }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// DataSource Implementation
//---------------------------------------------------------------------------//
class MyDataSource : public DataTransferKit::DataSource< std::vector<MyNode>,
							 std::vector<MyHex>,
							 std::vector<double> >
{
  private:

    std::vector<MyNode> d_nodes;
    std::vector<MyHex> d_elements;
    std::vector<double> d_element_data;
    MPI_Comm d_comm;

    void createMesh()
    {
	// Make some nodes.
	d_nodes.push_back( MyNode(0.0, 0.0, 0.0, 0) );
	d_nodes.push_back( MyNode(1.0, 0.0, 0.0, 4) );
	d_nodes.push_back( MyNode(1.0, 1.0, 0.0, 9) );
	d_nodes.push_back( MyNode(0.0, 1.0, 0.0, 2) );
	d_nodes.push_back( MyNode(0.0, 0.0, 1.0, 3) );
	d_nodes.push_back( MyNode(1.0, 0.0, 1.0, 8) );
	d_nodes.push_back( MyNode(1.0, 1.0, 1.0, 1) );
	d_nodes.push_back( MyNode(0.0, 1.0, 1.0, 6) );
	d_nodes.push_back( MyNode(0.0, 0.0, 2.0, 12) );
	d_nodes.push_back( MyNode(1.0, 0.0, 2.0, 7) );
	d_nodes.push_back( MyNode(1.0, 1.0, 2.0, 13) );
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 5) );

	// Make 2 hexahedrons.
	d_elements.push_back( MyHex( 0, 4, 9, 2, 3, 8, 1, 6, 0 ) );
	d_elements.push_back( MyHex( 3, 8, 1, 6, 12, 7, 13, 5, 1 ) ); 

	// Add some data for the hexes.
	d_element_data.push_back( 1.5 );
	d_element_data.push_back( 3.5 );
    }

  public:

    MyDataSource()
    { 
	// Build the mesh.
	createMesh();

	// Get the raw MPI_Comm out of Teuchos.
	Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
	Teuchos::RCP< const Teuchos::MpiComm<int> > mpi_comm = 
	    Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm<int> >( comm );
	Teuchos::RCP< const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm = 
	    mpi_comm->getRawMpiComm();
	d_comm = (*opaque_comm)();
   }

    ~MyDataSource()
    { /* ... */ }

    const MPI_Comm& getSourceComm()
    {
	return d_comm;
    }

    bool isFieldSupported( const std::string &field_name )
    {
	bool return_val = false;
	if ( field_name == "MY_DATA_FIELD" )
	{
	    return_val = true;
	}
	return return_val;
    }

    const std::vector<MyNode>& getSourceMeshNodes()
    {
	return d_nodes;
    }

    const std::vector<MyHex>& getSourceMeshElements()
    {
	return d_elements;
    }

    const std::vector<double> evaluateFieldOnTargetNodes( 
	const std::string &field_name,
	const std::vector<MyHex::handle_type> &element_handles,
	const std::vector<MyNode::coordinate_type> &node_coordinates )
    {
	if ( field_name == "MY_DATA_FIELD" )
	{
	    return d_element_data;
	}
	else
	{
	    std::vector<double> empty_vec;
	    return empty_vec;
	}
    }
};

//---------------------------------------------------------------------------//
// Mesh creation from data source.
//---------------------------------------------------------------------------//
namespace DataTransferKit
{

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

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

// DataSource test.
TEUCHOS_UNIT_TEST( Mesh, mesh_test )
{
    using namespace DataTransferKit;

    // Create a DataSource
    Teuchos::RCP< DataSource< std::vector<MyNode>,
			      std::vector<MyHex>,
			      std::vector<double> > > data_source 
	= Teuchos::rcp( new MyDataSource() );

    // Create a mesh.
    Teuchos::RCP< Mesh<MyHex::handle_type> > mesh = 
	createMeshFromDataSource( data_source );

    // Get the moab interface.
    moab::ErrorCode error;
    Mesh<MyHex::handle_type>::RCP_Moab moab = mesh->getMoab();
    
    // Grab the elements.
    moab::Range mesh_elements = mesh->getElements();

    // Check the moab mesh element data.
    moab::Range::const_iterator element_iterator;
    MyHex::handle_type native_handle = 0;
    for ( element_iterator = mesh_elements.begin();
	  element_iterator != mesh_elements.end();
	  ++element_iterator, ++native_handle )
    {
	TEST_ASSERT( mesh->getNativeHandle( *element_iterator ) == 
		     native_handle );

	TEST_ASSERT( moab->type_from_handle( *element_iterator ) ==
		     moab::MBHEX );
    }

    // Check the moab mesh vertex data.
    moab::Range connectivity;
    error = moab->get_connectivity( mesh_elements, connectivity );
    TEST_ASSERT( moab::MB_SUCCESS == error );

    std::vector<double> vertex_coords( 3 * connectivity.size() );
    error = moab->get_coords( connectivity, &vertex_coords[0] );
    TEST_ASSERT( moab::MB_SUCCESS == error );

    std::vector<double>::const_iterator moab_coord_iterator = 
	vertex_coords.begin();
    typename NodeTraits<MyNode>::const_coordinate_iterator node_coord_iterator;
    std::vector<MyNode> source_nodes = data_source->getSourceMeshNodes();
    std::vector<MyNode>::const_iterator node_iterator;

    for ( node_iterator = source_nodes.begin();
	  node_iterator != source_nodes.end();
	  ++node_iterator )
    {
	for ( node_coord_iterator = 
		  NodeTraits<MyNode>::coordsBegin( *node_iterator );
	      node_coord_iterator != 
		  NodeTraits<MyNode>::coordsEnd( *node_iterator );
	      ++node_coord_iterator, ++moab_coord_iterator )
	{
	    TEST_ASSERT( *node_coord_iterator == *moab_coord_iterator );
	}
    }
}

//---------------------------------------------------------------------------//
// end tstMesh.cpp
//---------------------------------------------------------------------------//

