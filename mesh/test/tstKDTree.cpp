//---------------------------------------------------------------------------//
/*! 
 * \file tstKDTree.cpp
 * \author Stuart R. Slattery
 * \brief kD-Tree unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_KDTree.hpp>
#include <DTK_RendezvousMesh.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_Exception.hpp>

#include <mpi.h>

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

    typedef int    handle_type;
    
    MyMesh() 
    { /* ... */ }

    MyMesh( const Teuchos::Array<int>& node_handles,
	    const Teuchos::Array<double>& coords,
	    const Teuchos::Array<int>& hex_handles,
	    const Teuchos::Array<int>& hex_connectivity )
	: d_node_handles( node_handles )
	, d_coords( coords )
	, d_hex_handles( hex_handles )
	, d_hex_connectivity( hex_connectivity )
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
    

  private:

    Teuchos::Array<int> d_node_handles;
    Teuchos::Array<double> d_coords;
    Teuchos::Array<int> d_hex_handles;
    Teuchos::Array<int> d_hex_connectivity;
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
    typedef Teuchos::Array<int>::const_iterator const_node_iterator;
    typedef Teuchos::Array<double>::const_iterator const_coordinate_iterator;
    typedef Teuchos::Array<int>::const_iterator const_element_iterator;    
    typedef Teuchos::Array<int>::const_iterator const_connectivity_iterator;    

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
   
    return MyMesh( node_handles, coords, hex_handles, hex_connectivity );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

// kD-Tree test.
TEUCHOS_UNIT_TEST( KDTree, kd_tree_test )
{
    using namespace DataTransferKit;

    // Create a mesh.
    MyMesh my_mesh = buildMyMesh();
    Teuchos::RCP< RendezvousMesh<MyMesh::handle_type> > mesh = 
	createRendezvousMesh( my_mesh );
    mesh->getMoab()->write_mesh( "tree.vtk" );
    // Create a KDTree.
    KDTree<MyMesh::handle_type> kd_tree( mesh );

    // Build the tree.
    kd_tree.build();

    // Search the tree for some points we know we will find.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 0.5;
    point_0[1] = 0.45;
    point_0[2] = 0.98;
    Teuchos::Array<double> point_1(3); 
    point_1[0] = 0.2;
    point_1[1] = 0.9;
    point_1[2] = 1.32;
    TEST_ASSERT( kd_tree.findPoint( point_0 ) == 0 );
    TEST_ASSERT( kd_tree.findPoint( point_1 ) == 1 );

    // Test points we know we won't find. We expect a PointNotFound()
    // exception to be thrown here.
    Teuchos::Array<double> point_2(3);
    point_2[0] = 2.9;
    point_2[1] = -0.5;
    point_2[2] = 9.5;
    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.1;
    point_3[1] = 1.5;
    point_3[2] = -4.8;
    bool caught_exception = false;
    try 
    {
	kd_tree.findPoint( point_2 );
    }
    catch( PointNotFound& exception )
    {
	caught_exception = true;
    }
    TEST_ASSERT( caught_exception );

    caught_exception = false;
    try 
    {
	kd_tree.findPoint( point_3 );
    }
    catch( PointNotFound& exception )
    {
	caught_exception = true;
    }
    TEST_ASSERT( caught_exception );
}

//---------------------------------------------------------------------------//
// end tstKDTree.cpp
//---------------------------------------------------------------------------//
