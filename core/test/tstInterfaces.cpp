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
#include <DTK_MeshTraits.hpp>
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
// Mesh Implementation
//---------------------------------------------------------------------------//

class MyMesh
{
  public:

    typedef int    handle_type;
    typedef double coordinate_type;
    
    MyMesh() 
    { /* ... */ }

    MyMesh( const std::vector<int>& node_handles,
	    const std::vector<double>& coords,
	    const std::vector<int>& quad_handles,
	    const std::vector<int>& quads_connectivity )
	: d_node_handles( node_handles )
	, d_coords( coords )
	, d_quad_handles( quad_handles )
	, d_quads_connectivity( quads_connectivity )
    { /* ... */ }

    ~MyMesh()
    { /* ... */ }

    std::vector<int>::const_iterator nodesBegin() const
    { return d_node_handles.begin(); }

    std::vector<int>::const_iterator nodesEnd() const
    { return d_node_handles.end(); }

    std::vector<double>::const_iterator coordsBegin() const
    { return d_coords.begin(); }

    std::vector<double>::const_iterator coordsEnd() const
    { return d_coords.end(); }

    std::vector<int>::const_iterator quadsBegin() const
    { return d_quad_handles.begin(); }

    std::vector<int>::const_iterator quadsEnd() const
    { return d_quad_handles.end(); }

    std::vector<int>::const_iterator connectivityBegin() const
    { return d_quads_connectivity.begin(); }

    std::vector<int>::const_iterator connectivityEnd() const
    { return d_quads_connectivity.end(); }
    

  private:

    std::vector<int> d_node_handles;
    std::vector<double> d_coords;
    std::vector<int> d_quad_handles;
    std::vector<int> d_quads_connectivity;
};

//---------------------------------------------------------------------------//
// DTK Traits Specializations
//---------------------------------------------------------------------------//
namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Mesh traits specialization for MyMesh
template<>
struct MeshTraits<MyMesh>
{
    typedef MyMesh::handle_type handle_type;
    typedef MyMesh::coordinate_type coordinate_type;
    typedef std::vector<int>::const_iterator const_handle_iterator;
    typedef std::vector<double>::const_iterator const_coordinate_iterator;
    

    static inline const_handle_iterator nodesBegin( const MyMesh& mesh )
    { return mesh.nodesBegin(); }

    static inline const_handle_iterator nodesEnd( const MyMesh& mesh )
    { return mesh.nodesEnd(); }

    static inline std::size_t nodeDim( const MyMesh& mesh )
    { return 3; }

    static inline bool interleavedNodeCoords( const MyMesh& mesh )
    { return true; }

    static inline const_coordinate_iterator coordsBegin( const MyMesh& mesh )
    { return mesh.coordsBegin(); }

    static inline const_coordinate_iterator coordsEnd( const MyMesh& mesh )
    { return mesh.coordsEnd(); }


    static inline std::size_t elementType( const MyMesh& mesh )
    { return DTK_FACE; }

    static inline std::size_t elementTopology( const MyMesh& mesh )
    { return DTK_QUADRILATERAL; }

    static inline std::size_t nodesPerElement( const MyMesh& mesh )
    { return 4; }

    static inline const_handle_iterator elementsBegin( const MyMesh& mesh )
    { return mesh.quadsBegin(); }

    static inline const_handle_iterator elementsEnd( const MyMesh& mesh )
    { return mesh.quadsEnd(); }

    static inline const_handle_iterator connectivityBegin( const MyMesh& mesh )
    { return mesh.connectivityBegin(); }

    static inline const_handle_iterator connectivityEnd( const MyMesh& mesh )
    { return mesh.connectivityEnd(); }
};

//---------------------------------------------------------------------------//
// FieldTraits specialization for the data and coordinate field.
template<>
struct FieldTraits< std::vector<double> >
{
    typedef double value_type;
    typedef std::vector<double>::iterator iterator;
    typedef std::vector<double>::const_iterator const_iterator;
    
    static inline std::size_t size( const std::vector<double>& data_field )
    { return data_field.size(); }

    static inline iterator begin( std::vector<double>& data_field )
    { return data_field.begin(); }

    static inline const_iterator begin( const std::vector<double>& data_field )
    { return data_field.begin(); }

    static inline iterator end( std::vector<double>& data_field )
    { return data_field.end(); }

    static inline const_iterator end( const std::vector<double>& data_field )
    { return data_field.end(); }

    static inline bool empty( const std::vector<double>& data_field )
    { return data_field.empty(); }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// DataSource Implementation
//---------------------------------------------------------------------------//
class MyDataSource : 
    public DataTransferKit::DataSource< MyMesh, std::vector<double> >
{
  private:

    MyMesh d_mesh;
    MPI_Comm d_comm;
    std::vector<double> d_data;

    void createMesh()
    {
	// Make some nodes.
	std::vector<int> node_handles;
	std::vector<double> node_coords;
	for ( int i = 0; i < 4; ++i )
	{
	    node_handles.push_back( i );
	    node_coords.push_back( 0.0 );
	    node_coords.push_back( 1.0 );
	    node_coords.push_back( 2.0 );
	}

	// Make a quadrilateral.
	std::vector<int> quad_handles;
	std::vector<int> quad_connectivity;
	quad_handles.push_back( 8 );
	for ( int i = 0; i < 4; ++i )
	{
	    quad_connectivity.push_back( i );
	}

	// Make a mesh.
	d_mesh = MyMesh( node_handles, node_coords,
			 quad_handles, quad_connectivity );

	// Add some data for the nodes.
	d_data.push_back( 1.5 );
	d_data.push_back( 3.5 );
	d_data.push_back( 5.5 );
	d_data.push_back( 7.5 );
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

    const MyMesh& getSourceMesh()
    {
	return d_mesh;
    }

    const std::vector<double> evaluateFieldOnTargetNodes( 
	const std::string &field_name,
	const std::vector<MyMesh::handle_type> &element_handles,
	const std::vector<MyMesh::coordinate_type> &node_coordinates )
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

class MyDataTarget : public DataTransferKit::DataTarget< std::vector<double>,
							 std::vector<double> >
{
  private:

    std::vector<double> d_coords;
    std::vector<double> d_data;
    MPI_Comm d_comm;

    void setup()
    {
	// Make some coordinates.
	for ( int i = 0; i < 4; ++i )
	{
	    d_coords.push_back( 0.0 );
	    d_coords.push_back( 1.0 );
	    d_coords.push_back( 2.0 );
	}

	// Allocate some memory for data for the nodes.
	d_data.resize( 4 );
    }

  public:

    MyDataTarget()
    { 
	// Setup.
	setup();

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

    bool interleavedCoordinates()
    {
	return true;
    }

    const std::vector<double>& getTargetCoordinates()
    {
	return d_coords;
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
// Tests
//---------------------------------------------------------------------------//
// DataSource test.
TEUCHOS_UNIT_TEST( DataSource, data_source_test )
{
    using namespace DataTransferKit;
    
    // Create a DataSource
    Teuchos::RCP< DataSource< MyMesh, std::vector<double> > > data_source = 
	Teuchos::rcp( new MyDataSource() );

    // Get the raw communicator and wrap it in a Teuchos::Comm interface.
    Teuchos::RCP< Teuchos::OpaqueWrapper<MPI_Comm> > raw_comm = 
	Teuchos::opaqueWrapper( data_source->getSourceComm() );
    Teuchos::RCP< Teuchos::Comm<int> > comm = 
	Teuchos::rcp( new Teuchos::MpiComm<int>( raw_comm ) );
    TEST_ASSERT( comm->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( comm->getSize() == getDefaultComm<int>()->getSize() );

    // Check that my data field is supported.
    TEST_ASSERT( data_source->isFieldSupported( "MY_DATA_FIELD" ) );

    // Get the mesh.
    MyMesh source_mesh = data_source->getSourceMesh();
    typename MeshTraits<MyMesh>::const_handle_iterator handle_iterator;
    typename MeshTraits<MyMesh>::const_coordinate_iterator coord_iterator;
    
    // Check the nodes.
    TEST_ASSERT( MeshTraits<MyMesh>::nodeDim( source_mesh ) == 3 );
    TEST_ASSERT( MeshTraits<MyMesh>::interleavedNodeCoords( source_mesh ) );
    TEST_ASSERT( std::distance( MeshTraits<MyMesh>::nodesBegin( source_mesh ),
				MeshTraits<MyMesh>::nodesEnd( source_mesh ) )
		 == 4 );
    
    int node_index = 0;
    for ( handle_iterator = MeshTraits<MyMesh>::nodesBegin( source_mesh );
	  handle_iterator != MeshTraits<MyMesh>::nodesEnd( source_mesh );
	  ++handle_iterator, ++node_index )
    {
	TEST_ASSERT( *handle_iterator == node_index );
    }

    for ( coord_iterator = MeshTraits<MyMesh>::coordsBegin( source_mesh );
	  coord_iterator != MeshTraits<MyMesh>::coordsEnd( source_mesh ); )
    {
	double coord_val = 0.0;
	for ( int i = 0; i < 3; ++i, ++coord_iterator )
	{
	    TEST_ASSERT( *coord_iterator == coord_val );
	    coord_val += 1.0;
	}
    }

    // Check the elements.
    TEST_ASSERT( MeshTraits<MyMesh>::elementType( source_mesh ) == DTK_FACE );
    TEST_ASSERT( MeshTraits<MyMesh>::elementTopology( source_mesh ) == 
		 DTK_QUADRILATERAL );
    TEST_ASSERT( MeshTraits<MyMesh>::nodesPerElement( source_mesh ) == 4 );
    TEST_ASSERT( std::distance( MeshTraits<MyMesh>::elementsBegin( source_mesh ),
				MeshTraits<MyMesh>::elementsEnd( source_mesh ) )
		 == 1 );
    TEST_ASSERT( std::distance( 
		     MeshTraits<MyMesh>::connectivityBegin( source_mesh ),
		     MeshTraits<MyMesh>::connectivityEnd( source_mesh ) )
		 == 4 );

    for ( handle_iterator = MeshTraits<MyMesh>::elementsBegin( source_mesh );
	  handle_iterator != MeshTraits<MyMesh>::elementsEnd( source_mesh );
	  ++handle_iterator )
    {
	TEST_ASSERT( *handle_iterator == 8 );
    }

    int conn_index = 0;
    for ( handle_iterator = 
	      MeshTraits<MyMesh>::connectivityBegin( source_mesh );
	  handle_iterator != 
	      MeshTraits<MyMesh>::connectivityEnd( source_mesh );
	  ++handle_iterator )
    {
	TEST_ASSERT( *handle_iterator == conn_index );
	++conn_index;
    }

    // Check the data.
    std::vector<MyMesh::handle_type> dummy_handles;
    std::vector<MyMesh::coordinate_type> dummy_coords;
    std::vector<double> data = 
	data_source->evaluateFieldOnTargetNodes( "MY_DATA_FIELD",
						 dummy_handles,
						 dummy_coords );

    TEST_ASSERT( FieldTraits< std::vector<double> >::size( data ) == 4 );
    double gold_data = 1.5;
    typename FieldTraits< std::vector<double> >::const_iterator data_iterator;
    for ( data_iterator = FieldTraits< std::vector<double> >::begin( data );
	  data_iterator != FieldTraits< std::vector<double> >::end( data );
	  ++data_iterator )
    {
	TEST_ASSERT( *data_iterator == gold_data );
	gold_data += 2.0;
    }
}

//---------------------------------------------------------------------------//
// DataTarget test.
TEUCHOS_UNIT_TEST( DataTarget, data_target_test )
{
    using namespace DataTransferKit;
    
    // Create a DataTarget
    Teuchos::RCP< DataTarget< std::vector<double>,
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

    // Check the target coordinates.
    TEST_ASSERT( data_target->interleavedCoordinates() );
    std::vector<double> target_coords = data_target->getTargetCoordinates();
    typename FieldTraits< std::vector<double> >::const_iterator coord_iterator;
    for ( coord_iterator = 
	      FieldTraits< std::vector<double> >::begin( target_coords );
	  coord_iterator != 
	      FieldTraits< std::vector<double> >::end( target_coords ); )
    {
	double coord_val = 0.0;
	for ( int i = 0; i < 3; ++i, ++coord_iterator )
	{
	    TEST_ASSERT( *coord_iterator == coord_val );
	    coord_val += 1.0;
	}
    }

    // Check that we can write mesh node data.
    std::vector<double> data = data_target->getTargetDataSpace( "MY_DATA_FIELD" );
    FieldTraits< std::vector<double> >::iterator write_iterator;
    double gold_data = 1.5;
    for ( write_iterator = FieldTraits< std::vector<double> >::begin( data );
	  write_iterator != FieldTraits< std::vector<double> >::end( data );
	  ++write_iterator )
    {
	*write_iterator = gold_data;
	gold_data += 1.0;
    }

    Teuchos::RCP<MyDataTarget> my_target = 
	Teuchos::rcp_dynamic_cast<MyDataTarget>( data_target );
    std::vector<double> target_data = my_target->getData();
    FieldTraits< std::vector<double> >::iterator read_iterator;
    gold_data = 1.5;
    for ( read_iterator = 
	      FieldTraits< std::vector<double> >::begin( target_data );
	  read_iterator != 
	      FieldTraits< std::vector<double> >::end( target_data );
	  ++read_iterator )
    {
	*read_iterator = gold_data;
	gold_data += 1.0;
    }
}

//---------------------------------------------------------------------------//
// Data copy test.
TEUCHOS_UNIT_TEST( DataSource, copy_test )
{
    using namespace DataTransferKit;

    // Create a DataSource
    Teuchos::RCP< DataSource< MyMesh, std::vector<double> > > data_source = 
	Teuchos::rcp( new MyDataSource() );

    // Create a DataTarget
    Teuchos::RCP< DataTarget< std::vector<double>,
			      std::vector<double> > > data_target = 
	Teuchos::rcp( new MyDataTarget() );

    // Copy from the source to the target.
    std::vector<int> dummy_handles;
    std::vector<double> dummy_coords;
    copyData( data_source->evaluateFieldOnTargetNodes( "MY_DATA_FIELD", 
						       dummy_handles, 
						       dummy_coords ), 
	      data_target->getTargetDataSpace( "MY_DATA_FIELD" ) );

    // Check the copy.
    Teuchos::RCP<MyDataTarget> my_target = 
	Teuchos::rcp_dynamic_cast<MyDataTarget>( data_target );
    std::vector<double> target_data = my_target->getData();
    FieldTraits< std::vector<double> >::iterator read_iterator;
    double gold_data = 1.5;
    for ( read_iterator = 
	      FieldTraits< std::vector<double> >::begin( target_data );
	  read_iterator != 
	      FieldTraits< std::vector<double> >::end( target_data );
	  ++read_iterator )
    {
	*read_iterator = gold_data;
	gold_data += 1.0;
    }
}

//---------------------------------------------------------------------------//
//                        end of tstInterfaces.cpp
//---------------------------------------------------------------------------//
