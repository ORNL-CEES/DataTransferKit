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

#include <DTK_DataSource.hpp>
#include <DTK_DataTarget.hpp>
#include <DTK_CoreTypes.hpp>
#include <DTK_NodeTraits.hpp>
#include <DTK_ElementTraits.hpp>
#include <DTK_FieldTraits.hpp>

#include <mpi.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

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

    int dim() const
    { return d_coords.size(); }

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

    typedef int handle_type;

    MyQuad( int node_1, int node_2, int node_3, int node_4, 
	    int handle )
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

// NodeTraits specialization for the MyNode implementation.
template<>
struct NodeTraits<MyNode>
{
    typedef typename MyNode::handle_type     handle_type;
    typedef typename MyNode::coordinate_type coordinate_type;
    typedef typename std::vector<double>::const_iterator coordinate_iterator;
    
    static inline std::size_t dim( const MyNode& node ) 
    { return node.dim();}
    
    static inline handle_type handle( const MyNode& node ) 
    { return node.handle(); }
    
    static inline coordinate_iterator coordsBegin( const MyNode& node ) 
    { return node.coordsBegin(); }

    static inline coordinate_iterator coordsEnd( const MyNode& node ) 
    { return node.coordsEnd(); }
};

// ElementTraits specialization for the MyQuad implementation.
template<>
struct ElementTraits<MyQuad>
{
    typedef typename MyQuad::handle_type handle_type;
    typedef typename std::vector<int>::const_iterator connectivity_iterator;

    static inline std::size_t topology()
    { return DTK_QUADRILATERAL; }

    static inline std::size_t num_nodes()
    { return 4; }

    static inline handle_type handle( const MyQuad &quad )
    { return quad.handle(); }

    static inline connectivity_iterator connectivityBegin( const MyQuad &quad )
    { return quad.connectivityBegin(); }

    static inline connectivity_iterator connectivityEnd( const MyQuad &quad )
    { return quad.connectivityEnd(); }
};

// FieldTraits specialization for the node field.
template<>
struct FieldTraits< std::vector<MyNode> >
{
    typedef MyNode value_type;
    typedef std::vector<MyNode>::const_iterator iterator;
    
    static inline std::size_t size( const std::vector<MyNode> &node_field )
    { return node_field.size(); }

    static iterator begin( const std::vector<MyNode> &node_field )
    { return node_field.begin(); }

    static inline iterator end( const std::vector<MyNode> &node_field )
    { return node_field.end(); }
};

// FieldTraits specialization for the element field.
template<>
struct FieldTraits< std::vector<MyQuad> >
{
    typedef MyQuad value_type;
    typedef std::vector<MyQuad>::const_iterator iterator;
    
    static inline std::size_t size( const std::vector<MyQuad> &quad_field )
    { return quad_field.size(); }

    static inline iterator begin( const std::vector<MyQuad> &quad_field )
    { return quad_field.begin(); }

    static inline iterator end( const std::vector<MyQuad> &quad_field )
    { return quad_field.end(); }
};

// FieldTraits specialization for the data field.
template<>
struct FieldTraits< std::vector<double> >
{
    typedef MyNode value_type;
    typedef std::vector<double>::iterator iterator;
    
    static inline std::size_t size( const std::vector<double> &data_field )
    { return data_field.size(); }

    static inline iterator begin( const std::vector<double> &data_field )
    { return data_field.begin(); }

    static inline iterator end( const std::vector<double> &data_field )
    { return data_field.end(); }
};

}

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
    std::vector<double> d_data;

    void createMesh()
    {
	// Make some nodes.
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 0) );
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 1) );
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 2) );
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 3) );

	// Make a quadrilateral.
	d_elements.push_back( MyQuad( 0, 1, 2, 3, 8 ) );

	// Add some data for the nodes.
	d_data.push_back( 1.5 );
	d_data.push_back( 3.5 );
	d_data.push_back( 7.5 );
	d_data.push_back( 9.5 );
    }

  public:

    MyDataSource()
    { 
	createMesh();
    }

    ~MyDataSource()
    { /* ... */ }

    const MPI_Comm& getSourceComm()
    {
	return getDefaultComm<int>()->getRawMpiComm();
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

    const std::vector<double>&
    getSourceNodeData( const std::string& field_name )
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
	createMesh();
    }

    ~MyDataTarget()
    { /* ... */ }

    const MPI_Comm& getTargetComm()
    {
	return getDefaultComm<int>()->getRawMpiComm();
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

    const std::vector<MyNode>& getTargetMeshNodes()
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
template<typename DataField>
void copyData( const DataField &source_field, DataField &target_field )
{
    using namespace DataTransferKit;

    TEST_ASSERT( FieldTraits<DataField>::size( source_field ) ==
		 FieldTraits<DataField>::size( target_field ) );

    typename FieldTraits<DataField>::iterator source_begin = 
	FieldTraits<DataField>::begin( source_field );

    typename FieldTraits<DataField>::iterator source_end = 
	FieldTraits<DataField>::end( source_field );

    typename FieldTraits<DataField>::iterator target_begin = 
	FieldTraits<DataField>::begin( target_field );

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

    TEST_ASSERT( FieldTraits<NodeField>::size( node_field ) == 3 );

    int node_index = 0;
    double coord_val = 0.0;
    typename FieldTraits<NodeField>::iterator node_iterator;
    for ( node_iterator = FieldTraits<NodeField>::begin( node_field );
	  node_iterator != FieldTraits<NodeField>::end( node_field );
	  ++node_iterator )
    {
	TEST_ASSERT( NodeTraits<NodeType>::dim( *node_iterator ) == 3 );
	TEST_ASSERT( NodeTraits<NodeType>::handle( *node_iterator ) == node_index );

	coord_val = 0.0;
	typename NodeTraits<NodeType>::coordinate_iterator coord_iterator;
	for ( coord_iterator = NodeTraits<NodeType>::coordsBegin( *node_iterator );
	      coord_iterator != NodeTraits<NodeType>::coordsEnd( *node_iterator );
	      ++coord_iterator )
	{
	    TEST_ASSERT( *coord_iterator == coord_val );
	    coord_val += 1.0;
	}

	++node_index;
    }
}

// Check the mesh elements.
template<typename ElementField>
void checkElements( const ElementField &element_field )
{
    using namespace DataTransferKit;

    typedef typename FieldTraits<ElementField>::value_type ElementType;

    TEST_ASSERT( FieldTraits<ElementField>::size( element_field ) == 1 );

    typename FieldTraits<ElementField>::iterator first_element = 
	FieldTraits<ElementField>::begin( element_field );

    TEST_ASSERT( ElementTraits<ElementType>::handle( *first_element ) == 8 );

    int conn_index = 0;
    typename ElementTraits<ElementType>::connectivity_iterator conn_iterator;
    for ( conn_iterator = ElementTraits<ElementType>::connectivityBegin( *first_element);
	  conn_iterator != ElementTraits<ElementType>::connectivityEnd( *first_element);
	  ++conn_iterator )
    {
	TEST_ASSERT( *conn_iterator ==  conn_index );
	++conn_index;
    }
}

// Check the node data.
template<typename DataField>
void checkNodeData( const DataField &data_field )
{
    using namespace DataTransferKit;

    TEST_ASSERT( FieldTraits<DataField>::size( data_field ) == 4 );

    double gold_data = 1.5;
    typename FieldTraits<DataField>::iterator data_iterator;
    for ( data_iterator = FieldTraits<DataField>::begin( data_field );
	  data_iterator != FieldTraits<DataField>::end( data_field );
	  ++data_iterator )
    {
	TEST_ASSERT( *data_iterator == gold_data );
	gold_data += 2.0;
    }
}

// Check that we can write data to the target.
template<typename DataField>
void checkWriteData( DataField &data_field )
{
    using namespace DataTransferKit;
    
    TEST_ASSERT( FieldTraits<DataField>::size( data_field ) == 4 );

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

    // Get the communicator and wrap it in a Teuchos::Comm interface.
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
    checkNodeData( data_source->getSourceNodeData( "MY_DATA_FIELD" ) );
}

// DataTarget test.
TEUCHOS_UNIT_TEST( DataTarget, data_target_test )
{
    // Create a DataTarget
    Teuchos::RCP< DataTarget< std::vector<MyNode>,
			      std::vector<double> > > data_target = 
	Teuchos::rcp( new MyDataTarget() );

    // Get the communicator and wrap it in a Teuchos::Comm interface.
    Teuchos::RCP< Teuchos::OpaqueWrapper<MPI_Comm> > raw_comm = 
	Teuchos::opaqueWrapper( data_target->getTargetComm() );
    Teuchos::RCP< Teuchos::Comm<int> > comm = 
	Teuchos::rcp( new Teuchos::MpiComm<int>( raw_comm ) );
    TEST_ASSERT( comm->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( comm->getSize() == getDefaultComm<int>()->getSize() );

    // Check that my data field is supported.
    TEST_ASSERT( data_target->isFieldSupported( "MY_DATA_FIELD" ) );

    // Check the mesh nodes.
    checkNodes( data_target->getTargetMeshNodes() );

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
    copyData( data_source->getSourceNodeData( "MY_DATA_FIELD" ), 
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
