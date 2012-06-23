//---------------------------------------------------------------------------//
/*!
 * \file tstMeshTools.cpp
 * \author Stuart R. Slattery
 * \brief Mesh tools unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <cstdlib>

#include <DTK_MeshTools.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_BoundingBox.hpp>

#include <mpi.h>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>

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

    typedef long int    global_ordinal_type;
    
    MyMesh() 
    { /* ... */ }

    MyMesh( const Teuchos::Array<global_ordinal_type>& node_handles,
	    const Teuchos::Array<double>& coords,
	    const Teuchos::Array<global_ordinal_type>& element_handles,
	    const Teuchos::Array<global_ordinal_type>& element_connectivity )
	: d_node_handles( node_handles )
	, d_coords( coords )
	, d_element_handles( element_handles )
	, d_element_connectivity( element_connectivity )
    { /* ... */ }

    ~MyMesh()
    { /* ... */ }

    Teuchos::Array<global_ordinal_type>::const_iterator nodesBegin() const
    { return d_node_handles.begin(); }

    Teuchos::Array<global_ordinal_type>::const_iterator nodesEnd() const
    { return d_node_handles.end(); }

    Teuchos::Array<double>::const_iterator coordsBegin() const
    { return d_coords.begin(); }

    Teuchos::Array<double>::const_iterator coordsEnd() const
    { return d_coords.end(); }

    Teuchos::Array<global_ordinal_type>::const_iterator elementsBegin() const
    { return d_element_handles.begin(); }

    Teuchos::Array<global_ordinal_type>::const_iterator elementsEnd() const
    { return d_element_handles.end(); }

    Teuchos::Array<global_ordinal_type>::const_iterator 
    connectivityBegin() const
    { return d_element_connectivity.begin(); }

    Teuchos::Array<global_ordinal_type>::const_iterator 
    connectivityEnd() const
    { return d_element_connectivity.end(); }
    

  private:

    Teuchos::Array<global_ordinal_type> d_node_handles;
    Teuchos::Array<double> d_coords;
    Teuchos::Array<global_ordinal_type> d_element_handles;
    Teuchos::Array<global_ordinal_type> d_element_connectivity;
};


//---------------------------------------------------------------------------//
// DTK implementations.
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
    typedef Teuchos::Array<global_ordinal_type>::const_iterator 
    const_node_iterator;
    typedef Teuchos::Array<double>::const_iterator 
    const_coordinate_iterator;
    typedef Teuchos::Array<global_ordinal_type>::const_iterator 
    const_element_iterator;
    typedef Teuchos::Array<global_ordinal_type>::const_iterator 
    const_connectivity_iterator;

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
    { return mesh.elementsBegin(); }

    static inline const_element_iterator elementsEnd( const MyMesh& mesh )
    { return mesh.elementsEnd(); }

    static inline const_connectivity_iterator connectivityBegin( const MyMesh& mesh )
    { return mesh.connectivityBegin(); }

    static inline const_connectivity_iterator connectivityEnd( const MyMesh& mesh )
    { return mesh.connectivityEnd(); }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Mesh create function.
//---------------------------------------------------------------------------//
MyMesh buildMyMesh( int my_rank, int my_size, int edge_length )
{
    // Make some nodes.
    int num_nodes = edge_length*edge_length*2;
    int node_dim = 3;
    Teuchos::Array<long int> node_handles( num_nodes );
    Teuchos::Array<double> coords( node_dim*num_nodes );
    int idx;
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length;
	    node_handles[ idx ] = (long int) num_nodes*my_rank + idx;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_nodes + idx ] = j;
	    coords[ 2*num_nodes + idx ] = 0.0;
	}
    }
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length + num_nodes / 2;
	    node_handles[ idx ] = (long int) num_nodes*my_rank + idx;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_nodes + idx ] = j;
	    coords[ 2*num_nodes + idx ] = 1.0;
	}
    }
    
    // Make the hexahedrons. 
    int num_elements = (edge_length-1)*(edge_length-1);
    Teuchos::Array<long int> hex_handles( num_elements );
    Teuchos::Array<long int> hex_connectivity( 8*num_elements );
    int elem_idx, node_idx;
    for ( int j = 0; j < (edge_length-1); ++j )
    {
	for ( int i = 0; i < (edge_length-1); ++i )
	{
	    node_idx = i + j*edge_length;
	    elem_idx = i + j*(edge_length-1);

	    hex_handles[elem_idx] = num_elements*my_rank + elem_idx;

	    hex_connectivity[elem_idx] 
		= node_handles[node_idx];

	    hex_connectivity[num_elements+elem_idx] 
		= node_handles[node_idx+1];

	    hex_connectivity[2*num_elements+elem_idx] 
		= node_handles[node_idx+edge_length+1];

	    hex_connectivity[3*num_elements+elem_idx] 
		= node_handles[node_idx+edge_length];

	    hex_connectivity[4*num_elements+elem_idx] 
		= node_handles[node_idx+num_nodes/2];

	    hex_connectivity[5*num_elements+elem_idx] 
		= node_handles[node_idx+num_nodes/2+1];

 	    hex_connectivity[6*num_elements+elem_idx] 
		= node_handles[node_idx+num_nodes/2+edge_length+1];

	    hex_connectivity[7*num_elements+elem_idx] 
		= node_handles[node_idx+num_nodes/2+edge_length];
	}
    }
    return MyMesh( node_handles, coords, hex_handles, hex_connectivity );
}

//---------------------------------------------------------------------------//
// Unit test.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( MeshTools, mesh_tools_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Setup source mesh.
    int edge_size = 4;
    MyMesh my_mesh = buildMyMesh( my_rank, my_size, edge_size );

    // Test the bounding boxes.
    BoundingBox local_box = MeshTools<MyMesh>::localBoundingBox( my_mesh );
    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == my_rank*(edge_size-1) );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == 0.0 );
    TEST_ASSERT( local_bounds[3] == (my_rank+1)*(edge_size-1) );
    TEST_ASSERT( local_bounds[4] == edge_size-1 );
    TEST_ASSERT( local_bounds[5] == 1.0 );

    BoundingBox global_box = 
	MeshTools<MyMesh>::globalBoundingBox( my_mesh, comm );
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == (my_size)*(edge_size-1) );
    TEST_ASSERT( global_bounds[4] == edge_size-1 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// end tstConsistentInterpolation2.cpp
//---------------------------------------------------------------------------//

