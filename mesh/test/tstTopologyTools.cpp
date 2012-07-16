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
#include <DTK_RendezvousMesh.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>

#include <MBInterface.hpp>
#include <MBCore.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
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

    typedef int    global_ordinal_type;
    
    MyMesh() 
    { /* ... */ }

    MyMesh( const Teuchos::Array<int>& node_handles,
	    const Teuchos::Array<double>& coords,
	    const Teuchos::Array<int>& hex_handles,
	    const Teuchos::Array<int>& hex_connectivity,
	    const Teuchos::Array<std::size_t>& permutation_list )
	: d_node_handles( node_handles )
	, d_coords( coords )
	, d_hex_handles( hex_handles )
	, d_hex_connectivity( hex_connectivity )
	, d_permutation_list( permutation_list )
    { /* ... */ }

    ~MyMesh()
    { /* ... */ }

    Teuchos::Array<int>::const_iterator nodesBegin() const
    { return d_node_handles.begin(); }

    Teuchos::Array<int>::const_iterator nodesEnd() const
    { return d_node_handles.end(); }

    Teuchos::Array<double>::const_iterator coordsBegin() const
    { return d_coords.begin(); }

    Teuchos::Array<double>::const_iterator coordsEnd() const
    { return d_coords.end(); }

    Teuchos::Array<int>::const_iterator hexesBegin() const
    { return d_hex_handles.begin(); }

    Teuchos::Array<int>::const_iterator hexesEnd() const
    { return d_hex_handles.end(); }

    Teuchos::Array<int>::const_iterator connectivityBegin() const
    { return d_hex_connectivity.begin(); }

    Teuchos::Array<int>::const_iterator connectivityEnd() const
    { return d_hex_connectivity.end(); }
    
    Teuchos::Array<std::size_t>::const_iterator permutationBegin() const
    { return d_permutation_list.begin(); }

    Teuchos::Array<std::size_t>::const_iterator permutationEnd() const
    { return d_permutation_list.end(); }


  private:

    Teuchos::Array<int> d_node_handles;
    Teuchos::Array<double> d_coords;
    Teuchos::Array<int> d_hex_handles;
    Teuchos::Array<int> d_hex_connectivity;
    Teuchos::Array<std::size_t> d_permutation_list;
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

    typedef MyMesh::global_ordinal_type global_ordinal_type;
    typedef Teuchos::Array<int>::const_iterator const_node_iterator;
    typedef Teuchos::Array<double>::const_iterator const_coordinate_iterator;
    typedef Teuchos::Array<int>::const_iterator const_element_iterator;
    typedef Teuchos::Array<int>::const_iterator const_connectivity_iterator;
    typedef Teuchos::Array<std::size_t>::const_iterator 
    const_permutation_iterator;


    static inline std::size_t nodeDim( const MyMesh& mesh )
    { return 3; }

    static inline const_node_iterator nodesBegin( const MyMesh& mesh )
    { return mesh.nodesBegin(); }

    static inline const_node_iterator nodesEnd( const MyMesh& mesh )
    { return mesh.nodesEnd(); }

    static inline const_coordinate_iterator coordsBegin( const MyMesh& mesh )
    { return mesh.coordsBegin(); }

    static inline const_coordinate_iterator coordsEnd( const MyMesh& mesh )
    { return mesh.coordsEnd(); }


    static inline std::size_t elementTopology( const MyMesh& mesh )
    { return DTK_HEXAHEDRON; }

    static inline std::size_t nodesPerElement( const MyMesh& mesh )
    { return 8; }

    static inline const_element_iterator elementsBegin( const MyMesh& mesh )
    { return mesh.hexesBegin(); }

    static inline const_element_iterator elementsEnd( const MyMesh& mesh )
    { return mesh.hexesEnd(); }

    static inline const_connectivity_iterator 
    connectivityBegin( const MyMesh& mesh )
    { return mesh.connectivityBegin(); }

    static inline const_connectivity_iterator 
    connectivityEnd( const MyMesh& mesh )
    { return mesh.connectivityEnd(); }

    static inline const_permutation_iterator permutationBegin( const MyMesh& mesh )
    { return mesh.permutationBegin(); }

    static inline const_permutation_iterator permutationEnd( const MyMesh& mesh )
    { return mesh.permutationEnd(); }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Mesh create funciton.
//---------------------------------------------------------------------------//
MyMesh buildMyMesh()
{
    // Make some nodes.
    Teuchos::Array<int> node_handles;
    Teuchos::Array<double> coords;

    // handles
    node_handles.push_back( 0 );
    node_handles.push_back( 4 );
    node_handles.push_back( 9 );
    node_handles.push_back( 2 );
    node_handles.push_back( 3 );
    node_handles.push_back( 8 );
    node_handles.push_back( 1 );
    node_handles.push_back( 6 );
    node_handles.push_back( 12 );
    node_handles.push_back( 7 );
    node_handles.push_back( 13 );
    node_handles.push_back( 5 );

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 );

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );

    // z
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );
    coords.push_back( 2.0 );
    coords.push_back( 2.0 );
    coords.push_back( 2.0 );
    coords.push_back( 2.0 );

    // Make 2 hexahedrons.
    Teuchos::Array<int> hex_handles;
    Teuchos::Array<int> hex_connectivity;
    
    // handles
    hex_handles.push_back( 0 );
    hex_handles.push_back( 1 );

    // 0
    hex_connectivity.push_back( 0 );
    hex_connectivity.push_back( 3 ); 

    // 1
    hex_connectivity.push_back( 4 ); 
    hex_connectivity.push_back( 8 );  

    // 2
    hex_connectivity.push_back( 9 );
    hex_connectivity.push_back( 1 ); 

    // 3
    hex_connectivity.push_back( 2 ); 
    hex_connectivity.push_back( 6 ); 

    // 4
    hex_connectivity.push_back( 3 );
    hex_connectivity.push_back( 12 ); 
   
    // 5
    hex_connectivity.push_back( 8 ); 
    hex_connectivity.push_back( 7 ); 

    // 6
    hex_connectivity.push_back( 1 ); 
    hex_connectivity.push_back( 13 ); 

    // 7
    hex_connectivity.push_back( 6 ); 
    hex_connectivity.push_back( 5 ); 
   
    Teuchos::Array<std::size_t> permutation_list( 8 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return MyMesh( node_handles, coords, hex_handles, hex_connectivity,
		   permutation_list );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Line point inclusion test.
TEUCHOS_UNIT_TEST( TopologyTools, line_test )
{
    using namespace DataTransferKit;

    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );
    error = moab->set_dimension( 1 );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    double vertex_0[1] = { -1.43 };
    double vertex_1[1] = { 2.98 };

    std::vector<moab::EntityHandle> vertices(2);
    error = moab->create_vertex( vertex_0, vertices[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_1, vertices[1] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    moab::EntityHandle line_segment;
    moab->create_element( moab::MBEDGE, &vertices[0], 2, line_segment );

    // Test the linear nodes function.
    TEST_ASSERT( TopologyTools::numLinearNodes( 
    		     moab->type_from_handle( line_segment ) ) == 2 );

    // Test the point inclusion test.
    Teuchos::Array<double> point_0(1);
    point_0[0] = 0.5;
    Teuchos::Array<double> point_1(1); 
    point_1[0] = 5.2;
    Teuchos::Array<double> point_2(1);
    point_2[0] = 2.98001;
    Teuchos::Array<double> point_3(1);
    point_3[0] = -1.5;
    Teuchos::Array<double> point_4(1);
    point_4[0] = 0.25;
    Teuchos::Array<double> point_5(1);
    point_5[0] = 2.98;

    TEST_ASSERT( TopologyTools::pointInElement( point_0, line_segment, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_1, line_segment, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, line_segment, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, line_segment, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_4, line_segment, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_5, line_segment, moab ) );
}

//---------------------------------------------------------------------------//
// Quad point inclusion test.
TEUCHOS_UNIT_TEST( TopologyTools, quad_test )
{
    using namespace DataTransferKit;

    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );
    error = moab->set_dimension( 2 );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    double vertex_0[2] = { -1.43, -2.3 };
    double vertex_1[2] = { 2.98, -2.12 };
    double vertex_2[2] = { 0.43, 4.2 };
    double vertex_3[2] = { -3.3, 1.1 };

    std::vector<moab::EntityHandle> vertices(3);
    error = moab->create_vertex( vertex_0, vertices[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_1, vertices[1] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_2, vertices[2] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_3, vertices[3] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    moab::EntityHandle quadrilateral;
    moab->create_element( moab::MBQUAD, &vertices[0], 4, quadrilateral );

    // Test the linear nodes function.
    TEST_ASSERT( TopologyTools::numLinearNodes( 
		     moab->type_from_handle( quadrilateral ) ) == 4 );

    // Test the point inclusion test.
    Teuchos::Array<double> point_0(2);
    point_0[0] = 0.5;
    point_0[1] = 0.45;
    Teuchos::Array<double> point_1(2); 
    point_1[0] = 5.2;
    point_1[1] = 0.9;
    Teuchos::Array<double> point_2(2);
    point_2[0] = 2.9;
    point_2[1] = -0.5;
    Teuchos::Array<double> point_3(2);
    point_3[0] = -11.5;
    point_3[1] = 1.5;
    Teuchos::Array<double> point_4(2);
    point_4[0] = 0.25;
    point_4[1] = -0.11;
    Teuchos::Array<double> point_5(2);
    point_5[0] = 0.43;
    point_5[1] = 4.2;

    TEST_ASSERT( TopologyTools::pointInElement( point_0, quadrilateral, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_1, quadrilateral, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, quadrilateral, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, quadrilateral, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_4, quadrilateral, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_5, quadrilateral, moab ) );
}

//---------------------------------------------------------------------------//
// Tri point inclusion test.
TEUCHOS_UNIT_TEST( TopologyTools, tri_test )
{
    using namespace DataTransferKit;

    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );
    error = moab->set_dimension( 2 );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    double vertex_0[2] = { -1.43, -2.3 };
    double vertex_1[2] = { 2.98, -2.12 };
    double vertex_2[2] = { 0.43, 4.2 };

    std::vector<moab::EntityHandle> vertices(3);
    error = moab->create_vertex( vertex_0, vertices[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_1, vertices[1] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_2, vertices[2] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    moab::EntityHandle triangle;
    moab->create_element( moab::MBTRI, &vertices[0], 3, triangle );

    // Test the linear nodes function.
    TEST_ASSERT( TopologyTools::numLinearNodes( 
		     moab->type_from_handle( triangle ) ) == 3 );

    // Test the point inclusion test.
    Teuchos::Array<double> point_0(2);
    point_0[0] = 0.5;
    point_0[1] = 0.45;
    Teuchos::Array<double> point_1(2); 
    point_1[0] = 5.2;
    point_1[1] = 0.9;
    Teuchos::Array<double> point_2(2);
    point_2[0] = 2.9;
    point_2[1] = -0.5;
    Teuchos::Array<double> point_3(2);
    point_3[0] = -1.5;
    point_3[1] = 1.5;
    Teuchos::Array<double> point_4(2);
    point_4[0] = 0.25;
    point_4[1] = -0.11;
    Teuchos::Array<double> point_5(2);
    point_5[0] = 0.43;
    point_5[1] = 4.2;

    TEST_ASSERT( TopologyTools::pointInElement( point_0, triangle, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_1, triangle, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, triangle, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, triangle, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_4, triangle, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_5, triangle, moab ) );
}

//---------------------------------------------------------------------------//
// Tet point inclusion test.
TEUCHOS_UNIT_TEST( TopologyTools, tet_test )
{
    using namespace DataTransferKit;

    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );

    double vertex_0[3] = { -1.43, -2.3, 3.4 };
    double vertex_1[3] = { 2.98, -2.12, 4.3 };
    double vertex_2[3] = { 0.43, 4.2, 2.4 };
    double vertex_3[3] = { 0.98, 3.77, 9.8 };

    std::vector<moab::EntityHandle> vertices(4);
    error = moab->create_vertex( vertex_0, vertices[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_1, vertices[1] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_2, vertices[2] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_3, vertices[3] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    moab::EntityHandle tetrahedron;
    moab->create_element( moab::MBTET, &vertices[0], 4, tetrahedron );

    // Test the linear nodes function.
    TEST_ASSERT( TopologyTools::numLinearNodes( 
		     moab->type_from_handle( tetrahedron ) ) == 4 );

    // Test the point inclusion test.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 0.5;
    point_0[1] = 0.45;
    point_0[2] = 6.2;
    Teuchos::Array<double> point_1(3); 
    point_1[0] = 0.2;
    point_1[1] = 0.9;
    point_1[2] = 1.32;
    Teuchos::Array<double> point_2(3);
    point_2[0] = 2.9;
    point_2[1] = -0.5;
    point_2[2] = 9.5;
    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.1;
    point_3[1] = 1.5;
    point_3[2] = -4.8;
    Teuchos::Array<double> point_4(3);
    point_4[0] = 0.25;
    point_4[1] = -0.11;
    point_4[2] = 5.4;
    Teuchos::Array<double> point_5(3);
    point_5[0] = 0.98;
    point_5[1] = 3.77;
    point_5[2] = 9.8;

    TEST_ASSERT( TopologyTools::pointInElement( point_0, tetrahedron, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_1, tetrahedron, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, tetrahedron, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, tetrahedron, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_4, tetrahedron, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_5, tetrahedron, moab ) );
}

//---------------------------------------------------------------------------//
// Pyramid point inclusion test.
TEUCHOS_UNIT_TEST( TopologyTools, pyramid_test )
{
    using namespace DataTransferKit;

    moab::ErrorCode error;
    Teuchos::RCP<moab::Interface> moab = Teuchos::rcp( new moab::Core() );

    double vertex_0[3] = { -1.43, -2.3, 3.4 };
    double vertex_1[3] = { 2.98, -2.12, 4.3 };
    double vertex_2[3] = { 0.43, 4.2, 2.4 };
    double vertex_3[3] = { -3.3, 1.1, 1.1 };
    double vertex_4[3] = { 0.98, 3.77, 9.8 };

    std::vector<moab::EntityHandle> vertices(5);
    error = moab->create_vertex( vertex_0, vertices[0] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_1, vertices[1] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_2, vertices[2] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_3, vertices[3] );
    TEST_ASSERT( error == moab::MB_SUCCESS );
    error = moab->create_vertex( vertex_4, vertices[4] );
    TEST_ASSERT( error == moab::MB_SUCCESS );

    moab::EntityHandle pyramid;
    moab->create_element( moab::MBPYRAMID, &vertices[0], 5, pyramid );

    // Test the linear nodes function.
    TEST_ASSERT( TopologyTools::numLinearNodes( 
		     moab->type_from_handle( pyramid ) ) == 5 );

    // Test the point inclusion test.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 0.5;
    point_0[1] = 0.45;
    point_0[2] = 6.2;
    Teuchos::Array<double> point_1(3); 
    point_1[0] = 0.2;
    point_1[1] = 0.9;
    point_1[2] = 1.32;
    Teuchos::Array<double> point_2(3);
    point_2[0] = 2.9;
    point_2[1] = -0.5;
    point_2[2] = 9.5;
    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.1;
    point_3[1] = 1.5;
    point_3[2] = -4.8;
    Teuchos::Array<double> point_4(3);
    point_4[0] = 0.25;
    point_4[1] = -0.11;
    point_4[2] = 5.2;
    Teuchos::Array<double> point_5(3);
    point_5[0] = 0.98;
    point_5[1] = 3.77;
    point_5[2] = 9.8;

    TEST_ASSERT( TopologyTools::pointInElement( point_0, pyramid, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_1, pyramid, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, pyramid, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, pyramid, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_4, pyramid, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_5, pyramid, moab ) );
}

//---------------------------------------------------------------------------//
// Hex point inclusion test.
TEUCHOS_UNIT_TEST( TopologyTools, hex_test )
{
    using namespace DataTransferKit;

    // Create a mesh.
    typedef MeshTraits<MyMesh> MT;
    Teuchos::ArrayRCP<MyMesh> mesh_blocks( 1 );
    mesh_blocks[0] = buildMyMesh();

    // Create a mesh manager.
    MeshManager<MyMesh> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );

    // Create a rendezvous mesh.
    Teuchos::RCP< RendezvousMesh<MyMesh::global_ordinal_type> > mesh = 
	createRendezvousMesh( mesh_manager );

    // Get the moab interface.
    RendezvousMesh<MyMesh::global_ordinal_type>::RCP_Moab moab = 
	mesh->getMoab();
    
    // Grab the elements.
    moab::Range mesh_elements = mesh->getElements();
    moab::EntityHandle hex_1 = mesh_elements[0];
    moab::EntityHandle hex_2 = mesh_elements[1];

    // Test the linear nodes function.
    TEST_ASSERT( 
	TopologyTools::numLinearNodes( moab->type_from_handle( hex_1 ) ) == 8 );
    TEST_ASSERT( 
	TopologyTools::numLinearNodes( moab->type_from_handle( hex_2 ) ) == 8 );

    // Test the point inclusion test.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 0.5;
    point_0[1] = 0.45;
    point_0[2] = 0.98;
    Teuchos::Array<double> point_1(3); 
    point_1[0] = 0.2;
    point_1[1] = 0.9;
    point_1[2] = 1.32;
    Teuchos::Array<double> point_2(3);
    point_2[0] = 2.9;
    point_2[1] = -0.5;
    point_2[2] = 9.5;
    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.1;
    point_3[1] = 1.5;
    point_3[2] = -4.8;
    Teuchos::Array<double> point_4(3);
    point_4[0] = 1.0;
    point_4[1] = 1.0;
    point_4[2] = 1.0;
    Teuchos::Array<double> point_5(3);
    point_5[0] = 0.0;
    point_5[1] = 0.0;
    point_5[2] = 1.0;

    TEST_ASSERT( TopologyTools::pointInElement( point_0, hex_1, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_1, hex_1, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, hex_1, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, hex_1, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_4, hex_1, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_5, hex_1, moab ) );

    TEST_ASSERT( !TopologyTools::pointInElement( point_0, hex_2, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_1, hex_2, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_2, hex_2, moab ) );
    TEST_ASSERT( !TopologyTools::pointInElement( point_3, hex_2, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_4, hex_2, moab ) );
    TEST_ASSERT( TopologyTools::pointInElement( point_5, hex_2, moab ) );
}

//---------------------------------------------------------------------------//
// end tstTopologyTools.cpp
//---------------------------------------------------------------------------//

