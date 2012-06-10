//---------------------------------------------------------------------------//
/*!
 * \file tstDataSource.cpp
 * \author Stuart R. Slattery
 * \brief DataSource unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_DataSource.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>

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
    
    MyMesh() 
    { /* ... */ }

    MyMesh( const std::vector<int>& node_handles,
	    const std::vector<double>& coords,
	    const std::vector<int>& hex_handles,
	    const std::vector<int>& hex_connectivity )
	: d_node_handles( node_handles )
	, d_coords( coords )
	, d_hex_handles( hex_handles )
	, d_hex_connectivity( hex_connectivity )
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

    std::vector<int>::const_iterator hexesBegin() const
    { return d_hex_handles.begin(); }

    std::vector<int>::const_iterator hexesEnd() const
    { return d_hex_handles.end(); }

    std::vector<int>::const_iterator connectivityBegin() const
    { return d_hex_connectivity.begin(); }

    std::vector<int>::const_iterator connectivityEnd() const
    { return d_hex_connectivity.end(); }
    

  private:

    std::vector<int> d_node_handles;
    std::vector<double> d_coords;
    std::vector<int> d_hex_handles;
    std::vector<int> d_hex_connectivity;
};

//---------------------------------------------------------------------------//
// DTK Traits Specializations
//---------------------------------------------------------------------------//
namespace DataTransferKit
{

//---------------------------------------------------------------------------//
// Mesh traits specialization for MyMesh
template<>
class MeshTraits<MyMesh>
{
  public:

    typedef MyMesh::handle_type handle_type;
    typedef std::vector<int>::const_iterator const_node_iterator;
    typedef std::vector<double>::const_iterator const_coordinate_iterator;
    typedef std::vector<int>::const_iterator const_element_iterator;    
    typedef std::vector<int>::const_iterator const_connectivity_iterator;    

    static inline const_node_iterator nodesBegin( const MyMesh& mesh )
    { return mesh.nodesBegin(); }

    static inline const_node_iterator nodesEnd( const MyMesh& mesh )
    { return mesh.nodesEnd(); }

    static inline const_coordinate_iterator coordsBegin( const MyMesh& mesh )
    { return mesh.coordsBegin(); }

    static inline const_coordinate_iterator coordsEnd( const MyMesh& mesh )
    { return mesh.coordsEnd(); }


    static inline std::size_t elementType( const MyMesh& mesh )
    { return DTK_REGION; }

    static inline std::size_t elementTopology( const MyMesh& mesh )
    { return DTK_HEXAHEDRON; }

    static inline std::size_t nodesPerElement( const MyMesh& mesh )
    { return 8; }

    static inline const_element_iterator elementsBegin( const MyMesh& mesh )
    { return mesh.hexesBegin(); }

    static inline const_element_iterator elementsEnd( const MyMesh& mesh )
    { return mesh.hexesEnd(); }

    static inline const_element_iterator connectivityBegin( const MyMesh& mesh )
    { return mesh.connectivityBegin(); }

    static inline const_element_iterator connectivityEnd( const MyMesh& mesh )
    { return mesh.connectivityEnd(); }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Mesh create funciton.
//---------------------------------------------------------------------------//
MyMesh buildMyMesh()
{
    // Make some nodes.
    std::vector<int> node_handles;
    std::vector<double> coords;

    node_handles.push_back( 0 );
    coords.push_back( 0.0 ); coords.push_back( 0.0 ); coords.push_back( 0.0 );

    node_handles.push_back( 4 );
    coords.push_back( 1.0 ); coords.push_back( 0.0 ); coords.push_back( 0.0 );

    node_handles.push_back( 9 );
    coords.push_back( 1.0 ); coords.push_back( 1.0 ); coords.push_back( 0.0 );

    node_handles.push_back( 2 );
    coords.push_back( 0.0 ); coords.push_back( 1.0 ); coords.push_back( 0.0 );

    node_handles.push_back( 3 );
    coords.push_back( 0.0 ); coords.push_back( 0.0 ); coords.push_back( 1.0 );

    node_handles.push_back( 8 );
    coords.push_back( 1.0 ); coords.push_back( 0.0 ); coords.push_back( 1.0 );

    node_handles.push_back( 1 );
    coords.push_back( 1.0 ); coords.push_back( 1.0 ); coords.push_back( 1.0 );

    node_handles.push_back( 6 );
    coords.push_back( 0.0 ); coords.push_back( 1.0 ); coords.push_back( 1.0 );

    node_handles.push_back( 12 );
    coords.push_back( 0.0 ); coords.push_back( 0.0 ); coords.push_back( 2.0 );

    node_handles.push_back( 7 );
    coords.push_back( 1.0 ); coords.push_back( 0.0 ); coords.push_back( 2.0 );

    node_handles.push_back( 13 );
    coords.push_back( 1.0 ); coords.push_back( 1.0 ); coords.push_back( 2.0 );

    node_handles.push_back( 5 );
    coords.push_back( 0.0 ); coords.push_back( 1.0 ); coords.push_back( 2.0 );

    // Make 2 hexahedrons.
    std::vector<int> hex_handles;
    std::vector<int> hex_connectivity;
    
    hex_handles.push_back( 0 );
    hex_connectivity.push_back( 0 ); hex_connectivity.push_back( 4 ); 
    hex_connectivity.push_back( 9 ); hex_connectivity.push_back( 2 ); 
    hex_connectivity.push_back( 3 ); hex_connectivity.push_back( 8 ); 
    hex_connectivity.push_back( 1 ); hex_connectivity.push_back( 6 ); 

    hex_handles.push_back( 1 );
    hex_connectivity.push_back( 3 ); hex_connectivity.push_back( 8 ); 
    hex_connectivity.push_back( 1 ); hex_connectivity.push_back( 6 ); 
    hex_connectivity.push_back( 12 ); hex_connectivity.push_back( 7 ); 
    hex_connectivity.push_back( 13 ); hex_connectivity.push_back( 5 ); 

    return MyMesh( node_handles, coords, hex_handles, hex_connectivity );
}

//---------------------------------------------------------------------------//
// FieldEvaluator implementation.
//---------------------------------------------------------------------------//
class MyEvaluator : 
    public DataTransferKit::DataSource<MyMesh,std::vector<double> >::FieldEvaluator
{
  public:

    typedef MyMesh::handle_type handle_type;

    MyEvaluator() 
    { /* ... */ }

    ~MyEvaluator() 
    { /* ... */ }
    
    std::vector<double> evaluate( const std::vector<handle_type>& elements,
				  const std::vector<double>& coords )
    {
	return std::vector<double>( elements.size(), 1.0 );
    }
};

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( DataSource, data_source_test )
{
    using namespace DataTransferKit;

    // Create a mesh.
    MyMesh my_mesh = buildMyMesh();

    // Get the raw MPI communicator.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    Teuchos::RCP< const Teuchos::MpiComm<int> > mpi_comm = 
	Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm<int> >( comm );
    Teuchos::RCP< const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm = 
	mpi_comm->getRawMpiComm();
    MPI_Comm raw_comm = (*opaque_comm)();

    // Create a data source object.
    typedef DataSource<MyMesh,std::vector<double> > DataSourceType;
    DataSourceType data_source( my_mesh, raw_comm );

    // Create and register the field evaluator.
    Teuchos::RCP<MyEvaluator> my_evaluator = Teuchos::rcp( new MyEvaluator() );
    data_source.registerFieldEvaluator( "MyField", my_evaluator );
    TEST_ASSERT( data_source.getFieldId( "MyField" ) == 0 );

    typename DataSourceType::RCP_FieldEvaluator evaluator = 
	data_source.getFieldEvaluator( 0 );
    std::vector<int> elements( 3, 1 );
    std::vector<double> coords( 9, 1.0 );
    std::vector<double> evaluated_data = evaluator->evaluate( elements, coords );
    TEST_ASSERT( elements.size() == evaluated_data.size() );

    std::vector<double>::const_iterator data_iterator;
    for ( data_iterator = evaluated_data.begin();
	  data_iterator != evaluated_data.end();
	  ++data_iterator )
    {
	TEST_ASSERT( *data_iterator == 1.0 );
    }
}

//---------------------------------------------------------------------------//
// end tstDataSource.cpp
//---------------------------------------------------------------------------//


