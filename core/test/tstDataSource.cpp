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

#include <DTK_DataSource.hpp>
#include <DTK_DataTarget.hpp>
#include <DTK_CoreTypes.hpp>
#include <DTK_NodeTraits.hpp>
#include <DTK_ElementTraits.hpp>
#include <DTK_FieldTraits.hpp>

#include "mpi.h"

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

    typedef int                                  handle_type;
    typedef double                               coordinate_type;
    typedef std::vector<double>::const_iterator  coordinate_const_iterator;
    
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

    coordinate_const_iterator coordsBegin() const
    { return d_coords.begin(); }

    coordinate_const_iterator coordsEnd() const
    { return d_coords.end(); }
};

// NodeTraits specialization for the MyNode implementation.
namespace DataTransferKit
{

template<>
struct NodeTraits<MyNode>
{
    typedef typename MyNode::handle_type               handle_type;
    typedef typename MyNode::coordinate_type           coordinate_type;
    typedef typename MyNode::coordinate_const_iterator coordinate_const_iterator;
    
    static inline std::size_t dim( const MyNode& node ) 
    { return node.dim();}
    
    static inline handle_type handle( const MyNode& node ) 
    { return node.handle(); }
    
    static inline coordinate_const_iterator 
    coordsBegin( const MyNode& node ) 
    { return node.coordsBegin(); }

    static inline coordinate_const_iterator 
    coordsEnd( const MyNode& node ) 
    { return node.coordsEnd(); }
};

}

//---------------------------------------------------------------------------//
// Element Implementation
//---------------------------------------------------------------------------//

class MyQuad
{
  private:

    std::size_t d_handle;
    std::vector<int> d_connectivity;

  public:

    typedef int                                           handle_type;
    typedef std::vector<int>::connectivity_const_iterator connectivity_const_iterator;

    MyQuad( int node_1, int node_2, int node_3, int node_4, 
	    int handle )
	, d_handle( handle )
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

    connectivity_const_iterator connectivityBegin() const
    { return d_connectivity.begin(); }

    connectivity_const_iterator connectivityEnd() const
    { return d_connectivity.end(); }
};

// ElementTraits specialization for the MyQuad implementation.
namespace DataTransferKit
{

template<>
struct ElementTraits<MyQuad>
{
    typedef typename MyQuad::handle_type                 handle_type;
    typedef typename MyQuad::connectivity_const_iterator connectivity_const_iterator;

    static inline std::size_t topology()
    { return DTK_QUADRILATERAL; }

    static inline handle_type handle( const MyQuad &quad )
    { return quad.handle(); }

    static inline connectivity_const_iterator
    connectivityBegin( const MyQuad &quad )
    { return quad.connectivityBegin(); }

    static inline connectivity_const_iterator
    connectivityEnd( const MyQuad &quad )
    { return quad.connectivityEnd(); }
};

}

//---------------------------------------------------------------------------//
// DataSource Implementation
//---------------------------------------------------------------------------//

class MyDataSource : public DataTransferKit::DataSource
{
  private:

    std::vector<MyNode> d_nodes;
    std::vector<MyElements> d_elements;
    std::vector<double> d_data;

    void createMesh()
    {
	// Make some nodes.
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 0) );
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 1) );
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 2) );
	d_nodes.push_back( MyNode(0.0, 1.0, 2.0, 3) );

	// Make a quadrilateral.
	d_elements.push_back( MyQuad( 0, 1, 2, 3, 0 ) );

	// Add some data for the nodes.
	d_data.push_back( 1.5 );
	d_data.push_back( 4.332 );
	d_data.push_back( 2.332 );
	d_data.push_back( 5.98 );
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
	std::vector<double> return_vec;
	if ( field_name == "MY_DATA_FIELD" )
	{
	    return_vec = d_data;
	}
	return return_vec;
    }
};

// FieldTraits specializations for the Node, Element, and Data fields.
namespace DataTransferKit
{

template<>
struct FieldTraits< std::vector<MyNode> >
{
    typedef MyNode                                          value_type;
    typedef typename::std::vector<MyNode>::const_iterator   const_iterator;
    
    static inline std::size_t size( const std::vector<MyNode> &node_field )
    { return node_field.size(); }

    static inline const_iterator begin( const std::vector<MyNode> &node_field )
    { return node_field.begin(); }

    static inline const_iterator end( const std::vector<MyNode> &node_field )
    { return node_field.end(); }
};

template<>
struct FieldTraits< std::vector<MyQuad> >
{
    typedef MyQuad                                         value_type;
    typedef typename std::vector<MyQuad>::const_iterator   const_iterator;
    
    static inline std::size_t size( const std::vector<MyQuad> &quad_field )
    { return quad_field.size(); }

    static inline const_iterator begin( const std::vector<MyQuad> &quad_field )
    { return quad_field.begin(); }

    static inline const_iterator end( const std::vector<MyQuad> &quad_field )
    { return quad_field.end(); }
};

template<>
struct FieldTraits< std::vector<double> >
{
    typedef MyNode                                         value_type;
    typedef typename std::vector<double>::const_iterator   const_iterator;
    
    static inline std::size_t size( const std::vector<double> &data_field )
    { return data_field.size(); }

    static inline const_iterator begin( const std::vector<double> &data_field )
    { return data_field.begin(); }

    static inline const_iterator end( const std::vector<double> &data_field )
    { return data_field.end(); }
};

}

//---------------------------------------------------------------------------//
// Check functions.
//---------------------------------------------------------------------------//
template<class NodeField>
void checkNodes( const NodeField &node_field )
{
    typedef typename FieldTraits<NodeField>::value_type NodeType;

    TEST_ASSERT( FieldTraits<NodeField>::size( node_field ) == 3 );

    int node_index = 0;
    double coord_val = 0.0;
    typename FieldTraits<NodeField>::const_iterator node_iterator;
    for ( node_iterator = FieldTraits<NodeField>::begin( node_field );
	  node_iterator != FieldTraits<NodeField>::end( node_field );
	  ++node_iterator )
    {
	TEST_ASSERT( NodeTraits<NodeType>::dim( *node_iterator ) == 3 );
	TEST_ASSERT( NodeTraits<NodeType>::handle( *node_iterator ) == node_index );

	coord_val = 0.0;
	typename NodeTraits<NodeType>::coordinate_const_iterator coord_iterator;
	for ( coord_iterator = NodeTraits<NodeType>::coordsBegin( *node_iterator );
	      coord_iterator != NodeTraits<NodeType>::coordsEnd( *node_iterator );
	      ++coord_iterator }
	{
	    TEST_ASSERT( *coord_iterator == coord_val );
	    coord_val += 1.0;
	}

	++node_index;
    }
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

namespace DataTransferKit {

TEUCHOS_UNIT_TEST( DataSource, data_source_test )
{
    // Create a DataSource
    Teuchos::RCP<DataSource> data_source = Teuchos::rcp( new MyDataSource() );

    // Get the communicator and wrap it in a Teuchos::Comm interface.
    Teuchos::RCP< Teuchos::OpaqueWrapper<MPI_Comm> > raw_comm = 
	Teuchos::opaqueWrapper( data_source->getSourceComm() );
    Teuchos::RCP< Teuchos::Comm<int> > comm = 
	Teuchos::rcp( new Teuchos::MpiComm( raw_comm ) );

    // Check that my data field is supported.
    TEST_ASSERT( data_source->isFieldSupported( "MY_DATA_FIELD" ) );

    // Check the mesh nodes.
    checkNodes( data_source->getSourceMeshNodes() );

    // Check the mesh elements.
    

    // Check the mesh node data.
}

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
//                        end of tstInterfaces.cpp
//---------------------------------------------------------------------------//
