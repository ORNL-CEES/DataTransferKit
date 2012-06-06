//---------------------------------------------------------------------------//
/*!
 * \file tstTopologyTools.cpp
 * \author Stuart R. Slattery
 * \brief Topology tools unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_TopologyTools.hpp>
#include <DTK_Mesh.hpp>
#include <DTK_CoreTypes.hpp>
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
struct MeshTraits<MyMesh>
{
    typedef MyMesh::handle_type handle_type;
    typedef std::vector<int>::const_iterator const_handle_iterator;
    typedef std::vector<double>::const_iterator const_coordinate_iterator;
    

    static inline const_handle_iterator nodesBegin( const MyMesh& mesh )
    { return mesh.nodesBegin(); }

    static inline const_handle_iterator nodesEnd( const MyMesh& mesh )
    { return mesh.nodesEnd(); }

    static inline bool interleavedCoords( const MyMesh& mesh )
    { return true; }

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

    static inline const_handle_iterator elementsBegin( const MyMesh& mesh )
    { return mesh.hexesBegin(); }

    static inline const_handle_iterator elementsEnd( const MyMesh& mesh )
    { return mesh.hexesEnd(); }

    static inline const_handle_iterator connectivityBegin( const MyMesh& mesh )
    { return mesh.connectivityBegin(); }

    static inline const_handle_iterator connectivityEnd( const MyMesh& mesh )
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
// Tests
//---------------------------------------------------------------------------//

// Point inclusion test.
TEUCHOS_UNIT_TEST( TopologyTools, topology_tools_test )
{
    using namespace DataTransferKit;

    // Create a mesh.
    MyMesh my_mesh = buildMyMesh();
    Teuchos::RCP< Mesh<MyMesh::handle_type> > mesh = createMesh( my_mesh );

    // Get the moab interface.
    Mesh<MyMesh::handle_type>::RCP_Moab moab = mesh->getMoab();
    
    // Grab the elements.
    moab::Range mesh_elements = mesh->getElements();
    moab::EntityHandle hex_1 = mesh_elements[0];
    moab::EntityHandle hex_2 = mesh_elements[1];

    // Test the linear nodes function.
    TEST_ASSERT( TopologyTools::numLinearNodes( moab->type_from_handle( hex_1 ) )
		 == 8 );
    TEST_ASSERT( TopologyTools::numLinearNodes( moab->type_from_handle( hex_2 ) )
		 == 8 );

    // Test the point inclusion test.
    double point_0[3] = { 0.5, 0.45, 0.98 };
    double point_1[3] = { 0.2, 0.9, 1.32 };
    double point_2[3] = { 2.9, -0.5, 9.5 };
    double point_3[3] = { 0.1, 1.5, -4.8 };

    TEST_ASSERT( TopologyTools::pointInElement( point_0, hex_1, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_1, hex_1, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, hex_1, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, hex_1, moab ) );

    TEST_ASSERT( !TopologyTools::pointInElement( point_0, hex_2, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_1, hex_2, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, hex_2, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, hex_2, moab ) );
}

//---------------------------------------------------------------------------//
// end tstTopologyTools.cpp
//---------------------------------------------------------------------------//

