//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstInterfaces.cpp
 * \author Stuart Slattery
 * \date   Thu Dec 01 16:50:04 2011
 * \brief  Unit tests for the data transfer pure virtual interfaces.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_DataSource.hpp>
#include <DTK_DataTarget.hpp>
#include <DTK_CoreTypes.hpp>
#include <DTK_NodeTraits.hpp>
#include <DTK_ElementTraits.hpp>
#include <DTK_FieldTraits.hpp>

#include <mpi.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>

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

class MyQuad
{
  private:

    std::size_t d_handle;
    std::vector<int> d_connectivity;

  public:

    typedef std::size_t handle_type;

    MyQuad( int node_1, int node_2, int node_3, int node_4, 
	    std::size_t handle )
	: d_handle( handle )
    {
	d_connectivity.push_back( node_1 );
	d_connectivity.push_back( node_2 );
	d_connectivity.push_back( node_3 );
	d_connectivity.push_back( node_4 );
    }

    ~MyQuad()
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
// ElementTraits specialization for the MyQuad implementation.
template<>
struct ElementTraits<MyQuad>
{
    typedef typename MyQuad::handle_type              handle_type;
    typedef typename std::vector<int>::const_iterator 
    const_connectivity_iterator;

    static inline std::size_t type()
    { return DTK_FACE; }

    static inline std::size_t topology()
    { return DTK_QUADRILATERAL; }

    static inline std::size_t numNodes()
    { return 4; }

    static inline handle_type handle( const MyQuad &quad )
    { return quad.handle(); }

    static inline const_connectivity_iterator 
    connectivityBegin( const MyQuad &quad )
    { return quad.connectivityBegin(); }

    static inline const_connectivity_iterator 
    connectivityEnd( const MyQuad &quad )
    { return quad.connectivityEnd(); }
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
struct FieldTraits< std::vector<MyQuad> >
{
    typedef MyQuad                                value_type;
    typedef std::vector<MyQuad>::iterator         iterator;
    typedef std::vector<MyQuad>::const_iterator   const_iterator;
    
    static inline std::size_t size( const std::vector<MyQuad> &quad_field )
    { return quad_field.size(); }

    static inline iterator begin( std::vector<MyQuad> &quad_field )
    { return quad_field.begin(); }

    static inline const_iterator begin( const std::vector<MyQuad> &quad_field )
    { return quad_field.begin(); }

    static inline iterator end( std::vector<MyQuad> &quad_field )
    { return quad_field.end(); }

    static inline const_iterator end( const std::vector<MyQuad> &quad_field )
    { return quad_field.end(); }

    static inline bool empty(  const std::vector<MyQuad> &quad_field )
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
							 std::vector<MyQuad>,
							 std::vector<double> >
{
  private:

    std::vector<MyNode> d_nodes;
    std::vector<MyQuad> d_elements;
    std::vector<double> d_element_data;
    MPI_Comm d_comm;

    void createMesh()
    {
	// Make some nodes.
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 0) );
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 1) );
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 2) );
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 3) );

	// Make a quadrilateral.
	d_elements.push_back( MyQuad( 0, 1, 2, 3, 8 ) );

	// Add some data for the elements.
	d_element_data.push_back( 1.5 );
	d_element_data.push_back( 3.5 );
	d_element_data.push_back( 5.5 );
	d_element_data.push_back( 7.5 );
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

    const std::vector<MyQuad>& getSourceMeshElements()
    {
	return d_elements;
    }

    const std::vector<double> evaluateFieldOnTargetNodes( 
	const std::string &field_name,
	const std::vector<MyQuad::handle_type> &element_handles,
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
// DataTarget implementation.
//---------------------------------------------------------------------------//

class MyDataTarget : public DataTransferKit::DataTarget< std::vector<MyNode>,
							 std::vector<double> >
{
  private:

    std::vector<MyNode> d_nodes;
    std::vector<MyQuad> d_elements;
    std::vector<double> d_data;
    MPI_Comm d_comm;

    void createMesh()
    {
	// Make some nodes.
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 0) );
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 1) );
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 2) );
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 3) );

	// Make a quadrilateral.
	d_elements.push_back( MyQuad( 0, 1, 2, 3, 8 ) );

	// Allocate some memory for data for the nodes.
	d_data.resize( 4 );
    }

  public:

    MyDataTarget()
    { 
	// Build the mesh.
	createMesh();

	// Get the raw MPI_Comm out of Teuchos.
	Teuchos::RCP< const Teuchos::MpiComm<int> > mpi_comm = 
	    Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm<int> >(
		getDefaultComm<int>() );
	Teuchos::RCP< const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm = 
	    mpi_comm->getRawMpiComm();
	d_comm = (*opaque_comm)();
    }

    ~MyDataTarget()
    { /* ... */ }

    const MPI_Comm& getTargetComm()
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

    const std::vector<MyNode>& getTargetNodes()
    {
	return d_nodes;
    }

    std::vector<double>&
    getTargetDataSpace( const std::string& field_name )
    {
	if ( field_name == "MY_DATA_FIELD" )
	{
	    return d_data;
	}
	else
	{
	    std::vector<double> empty_vec;
	    return empty_vec;
	}
    }

    const std::vector<double>& getData() const
    { return d_data; }
};

//---------------------------------------------------------------------------//
// Copy function.
//---------------------------------------------------------------------------//
template<typename SourceDataField, typename TargetDataField>
void copyData( const SourceDataField &source_field, 
	       TargetDataField &target_field )
{
    using namespace DataTransferKit;

    typedef typename FieldTraits<SourceDataField>::value_type source_type;
    typedef typename FieldTraits<TargetDataField>::value_type target_type;

    bool same_type = 
	Teuchos::TypeTraits::is_same<source_type,target_type>::value;
    assert( same_type );

    assert( FieldTraits<SourceDataField>::size( source_field ) ==
		 FieldTraits<TargetDataField>::size( target_field ) );

    typename FieldTraits<SourceDataField>::const_iterator source_begin = 
	FieldTraits<SourceDataField>::begin( source_field );

    typename FieldTraits<SourceDataField>::const_iterator source_end = 
	FieldTraits<SourceDataField>::end( source_field );

    typename FieldTraits<TargetDataField>::iterator target_begin = 
	FieldTraits<TargetDataField>::begin( target_field );

    std::copy( source_begin, source_end, target_begin );
}

//---------------------------------------------------------------------------//
// Check functions.
//---------------------------------------------------------------------------//
// Check the mesh nodes.
template<typename NodeField>
void checkNodes( const NodeField &node_field )
{
    using namespace DataTransferKit;

    typedef typename FieldTraits<NodeField>::value_type NodeType;

    assert( FieldTraits<NodeField>::size( node_field ) == 4 );

    int node_index = 0;
    double coord_val = 0.0;
    typename FieldTraits<NodeField>::const_iterator node_iterator;
    for ( node_iterator = FieldTraits<NodeField>::begin( node_field );
	  node_iterator != FieldTraits<NodeField>::end( node_field );
	  ++node_iterator )
    {
	assert( NodeTraits<NodeType>::dim() == 3 );
	assert( NodeTraits<NodeType>::handle( *node_iterator ) == node_index );

	coord_val = 0.0;
	typename NodeTraits<NodeType>::const_coordinate_iterator coord_iterator;
	for ( coord_iterator = NodeTraits<NodeType>::coordsBegin( *node_iterator );
	      coord_iterator != NodeTraits<NodeType>::coordsEnd( *node_iterator );
	      ++coord_iterator )
	{
	    assert( *coord_iterator == coord_val );
	    coord_val += 1.0;
	}

	++node_index;
    }
}

//---------------------------------------------------------------------------//
// Check the mesh elements.
template<typename ElementField>
void checkElements( const ElementField &element_field )
{
    using namespace DataTransferKit;

    typedef typename FieldTraits<ElementField>::value_type ElementType;

    assert( FieldTraits<ElementField>::size( element_field ) == 1 );

    typename FieldTraits<ElementField>::const_iterator first_element = 
	FieldTraits<ElementField>::begin( element_field );

    assert( ElementTraits<ElementType>::handle( *first_element ) == 8 );

    int conn_index = 0;
    typename ElementTraits<ElementType>::const_connectivity_iterator conn_iterator;
    for ( conn_iterator = ElementTraits<ElementType>::connectivityBegin( *first_element);
	  conn_iterator != ElementTraits<ElementType>::connectivityEnd( *first_element);
	  ++conn_iterator )
    {
	assert( *conn_iterator ==  conn_index );
	++conn_index;
    }
}

//---------------------------------------------------------------------------//
// Check the node data.
template<typename DataField>
void checkNodeData( const DataField &data_field )
{
    using namespace DataTransferKit;

    assert( FieldTraits<DataField>::size( data_field ) == 4 );

    double gold_data = 1.5;
    typename FieldTraits<DataField>::const_iterator data_iterator;
    for ( data_iterator = FieldTraits<DataField>::begin( data_field );
	  data_iterator != FieldTraits<DataField>::end( data_field );
	  ++data_iterator )
    {
	assert( *data_iterator == gold_data );
	gold_data += 2.0;
    }
}

//---------------------------------------------------------------------------//
// Check that we can write data to the target.
template<typename DataField>
void checkWriteData( DataField &data_field )
{
    using namespace DataTransferKit;
    
    assert( FieldTraits<DataField>::size( data_field ) == 4 );

    double gold_data = 1.5;
    typename FieldTraits<DataField>::iterator data_iterator;
    for ( data_iterator = FieldTraits<DataField>::begin( data_field );
	  data_iterator != FieldTraits<DataField>::end( data_field );
	  ++data_iterator )
    {
	*data_iterator = gold_data;
	gold_data += 2.0;
    }
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
namespace DataTransferKit {

// DataSource test.
TEUCHOS_UNIT_TEST( DataSource, data_source_test )
{
    // Create a DataSource
    Teuchos::RCP< DataSource< std::vector<MyNode>,
			      std::vector<MyQuad>,
			      std::vector<double> > > data_source 
			      = Teuchos::rcp( new MyDataSource() );

    // Get the raw communicator and wrap it in a Teuchos::Comm interface.
    Teuchos::RCP< Teuchos::OpaqueWrapper<MPI_Comm> > raw_comm = 
	Teuchos::opaqueWrapper( data_source->getSourceComm() );
    Teuchos::RCP< Teuchos::Comm<int> > comm = 
	Teuchos::rcp( new Teuchos::MpiComm<int>( raw_comm ) );
    TEST_ASSERT( comm->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( comm->getSize() == getDefaultComm<int>()->getSize() );

    // Check that my data field is supported.
    TEST_ASSERT( data_source->isFieldSupported( "MY_DATA_FIELD" ) );

    // Check the mesh nodes.
    checkNodes( data_source->getSourceMeshNodes() );

    // Check the mesh elements.
    checkElements( data_source->getSourceMeshElements() );

    // Check the mesh node data.
    std::vector<MyQuad::handle_type> dummy_handles;
    std::vector<MyNode::coordinate_type> dummy_coords;
    checkNodeData( data_source->evaluateFieldOnTargetNodes( "MY_DATA_FIELD",
							    dummy_handles,
							    dummy_coords ) );
}

// DataTarget test.
TEUCHOS_UNIT_TEST( DataTarget, data_target_test )
{
    // Create a DataTarget
    Teuchos::RCP< DataTarget< std::vector<MyNode>,
			      std::vector<double> > > data_target = 
	Teuchos::rcp( new MyDataTarget() );

    // Get the raw communicator and wrap it in a Teuchos::Comm interface.
    Teuchos::RCP< Teuchos::OpaqueWrapper<MPI_Comm> > raw_comm = 
	Teuchos::opaqueWrapper( data_target->getTargetComm() );
    Teuchos::RCP< Teuchos::Comm<int> > comm = 
	Teuchos::rcp( new Teuchos::MpiComm<int>( raw_comm ) );
    TEST_ASSERT( comm->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( comm->getSize() == getDefaultComm<int>()->getSize() );

    // Check that my data field is supported.
    TEST_ASSERT( data_target->isFieldSupported( "MY_DATA_FIELD" ) );

    // Check the mesh nodes.
    checkNodes( data_target->getTargetNodes() );

    // Check that we can write mesh node data.
    checkWriteData( data_target->getTargetDataSpace( "MY_DATA_FIELD" ) );
    Teuchos::RCP<MyDataTarget> my_target = 
	Teuchos::rcp_dynamic_cast<MyDataTarget>( data_target );
    checkNodeData( my_target->getData() );
}

// Data copy test.
TEUCHOS_UNIT_TEST( DataSource, copy_test )
{
    // Create a DataSource
    Teuchos::RCP< DataSource< std::vector<MyNode>,
			      std::vector<MyQuad>,
			      std::vector<double> > > data_source = 
	Teuchos::rcp( new MyDataSource() );

    // Create a DataTarget
    Teuchos::RCP< DataTarget< std::vector<MyNode>,
			      std::vector<double> > > data_target = 
	Teuchos::rcp( new MyDataTarget() );

    // Copy from the source to the target.
    std::vector<MyQuad::handle_type> dummy_handles;
    std::vector<MyNode::coordinate_type> dummy_coords;
    copyData( data_source->evaluateFieldOnTargetNodes( "MY_DATA_FIELD", 
						       dummy_handles, 
						       dummy_coords ), 
	      data_target->getTargetDataSpace( "MY_DATA_FIELD" ) );

    // Check the copy.
    Teuchos::RCP<MyDataTarget> my_target = 
	Teuchos::rcp_dynamic_cast<MyDataTarget>( data_target );
    checkNodeData( my_target->getData() );
}

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
//                        end of tstInterfaces.cpp
//---------------------------------------------------------------------------//
