//---------------------------------------------------------------------------//
/*!
 * \file tstRendezvousMesh.cpp
 * \author Stuart R. Slattery
 * \brief RendezvousMesh unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_RendezvousMesh.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
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
	    const Teuchos::Array<int>& hex_connectivity )
	: d_node_handles( node_handles )
	, d_coords( coords )
	, d_hex_handles( hex_handles )
	, d_hex_connectivity( hex_connectivity )
	, d_permutation_list( 8 )
    {
	for ( int i = 0; i < d_permutation_list.size(); ++i )
	{
	    d_permutation_list[i] = i;
	}
    }

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

    static inline const_connectivity_iterator connectivityBegin( const MyMesh& mesh )
    { return mesh.connectivityBegin(); }

    static inline const_connectivity_iterator connectivityEnd( const MyMesh& mesh )
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
   
    return MyMesh( node_handles, coords, hex_handles, hex_connectivity );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( RendezvousMesh, rendezvous_mesh_test )
{
    using namespace DataTransferKit;

    // Create a mesh.
    typedef MeshTraits<MyMesh> MT;
    MyMesh my_mesh = buildMyMesh();
    Teuchos::RCP< RendezvousMesh<MyMesh::global_ordinal_type> > mesh = 
	createRendezvousMesh( my_mesh );

    // Get the moab interface.
    moab::ErrorCode error;
    RendezvousMesh<MyMesh::global_ordinal_type>::RCP_Moab moab = mesh->getMoab();
    
    // Grab the elements.
    moab::Range mesh_elements = mesh->getElements();

    // Check the moab mesh element data.
    moab::Range::const_iterator element_iterator;
    MyMesh::global_ordinal_type native_handle = 0;
    for ( element_iterator = mesh_elements.begin();
	  element_iterator != mesh_elements.end();
	  ++element_iterator, ++native_handle )
    {
	TEST_ASSERT( mesh->getNativeOrdinal( *element_iterator ) == 
		     native_handle );

	TEST_ASSERT( moab->type_from_handle( *element_iterator ) ==
		     moab::MBHEX );
    }

    // Check the moab mesh vertex data.
    moab::Range connectivity;
    error = moab->get_connectivity( mesh_elements, connectivity );
    TEST_ASSERT( moab::MB_SUCCESS == error );

    Teuchos::Array<double> vertex_coords( 3 * connectivity.size() );
    error = moab->get_coords( connectivity, &vertex_coords[0] );
    TEST_ASSERT( moab::MB_SUCCESS == error );

    int num_nodes = connectivity.size();
    Teuchos::Array<double>::const_iterator moab_coord_iterator;
    typename MT::const_coordinate_iterator coord_iterator = 
	MT::coordsBegin( my_mesh );
    int i = 0;
    for ( moab_coord_iterator = vertex_coords.begin();
	  moab_coord_iterator != vertex_coords.end(); ++i )
    {
	for ( int d = 0; d < 3; ++d, ++moab_coord_iterator )
	{
	    TEST_ASSERT( coord_iterator[d*num_nodes + i] == 
			 *moab_coord_iterator );
	}
    }
}

//---------------------------------------------------------------------------//
// end tstRendezvousMesh.cpp
//---------------------------------------------------------------------------//

