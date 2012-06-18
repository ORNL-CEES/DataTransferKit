//---------------------------------------------------------------------------//
/*!
 * \file tstRendezvous.cpp
 * \author Stuart R. Slattery
 * \brief Rendezvous unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_Rendezvous.hpp>
#include <DTK_RendezvousMesh.hpp>
#include <DTK_BoundingBox.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
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
	    const Teuchos::Array<int>& quad_handles,
	    const Teuchos::Array<int>& quad_connectivity )
	: d_node_handles( node_handles )
	, d_coords( coords )
	, d_quad_handles( quad_handles )
	, d_quad_connectivity( quad_connectivity )
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

    Teuchos::Array<int>::const_iterator quadsBegin() const
    { return d_quad_handles.begin(); }

    Teuchos::Array<int>::const_iterator quadsEnd() const
    { return d_quad_handles.end(); }

    Teuchos::Array<int>::const_iterator connectivityBegin() const
    { return d_quad_connectivity.begin(); }

    Teuchos::Array<int>::const_iterator connectivityEnd() const
    { return d_quad_connectivity.end(); }
    

  private:

    Teuchos::Array<int> d_node_handles;
    Teuchos::Array<double> d_coords;
    Teuchos::Array<int> d_quad_handles;
    Teuchos::Array<int> d_quad_connectivity;
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

    static inline std::size_t nodeDim( const MyMesh& mesh )
    { return 2; }

    static inline const_node_iterator nodesBegin( const MyMesh& mesh )
    { return mesh.nodesBegin(); }

    static inline const_node_iterator nodesEnd( const MyMesh& mesh )
    { return mesh.nodesEnd(); }

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


    static inline const_element_iterator elementsBegin( const MyMesh& mesh )
    { return mesh.quadsBegin(); }

    static inline const_element_iterator elementsEnd( const MyMesh& mesh )
    { return mesh.quadsEnd(); }

    static inline const_connectivity_iterator connectivityBegin( const MyMesh& mesh )
    { return mesh.connectivityBegin(); }

    static inline const_connectivity_iterator connectivityEnd( const MyMesh& mesh )
    { return mesh.connectivityEnd(); }
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Mesh create funciton.
//---------------------------------------------------------------------------//
MyMesh buildMyMesh()
{
    int my_rank = getDefaultComm<int>()->getRank();

    // Make some nodes.
    int num_nodes = 10;
    int node_dim = 2;
    Teuchos::Array<int> node_handles( num_nodes );
    Teuchos::Array<double> coords( node_dim*num_nodes );

    for ( int i = 0; i < num_nodes; ++i )
    {
	node_handles[i] = (num_nodes / 2)*my_rank + i;
    }
    for ( int i = 0; i < num_nodes / 2; ++i )
    {
	coords[ i ] = my_rank;
	coords[ num_nodes + i ] = i;
    }
    for ( int i = num_nodes / 2; i < num_nodes; ++i )
    {
	coords[ i ] = my_rank + 1;
	coords[ num_nodes + i ] = i - num_nodes/2;
    }
    
    // Make the quads.
    int num_quads = 4;
    Teuchos::Array<int> quad_handles( num_quads );
    Teuchos::Array<int> quad_connectivity( 4*num_quads );
    
    for ( int i = 0; i < num_quads; ++i )
    {
	quad_handles[ i ] = num_quads*my_rank + i;
	quad_connectivity[ i ] = node_handles[i];
	quad_connectivity[ num_quads + i ] = node_handles[num_nodes/2 + i];
	quad_connectivity[ 2*num_quads + i ] = node_handles[num_nodes/2 + i + 1];
	quad_connectivity[ 3*num_quads + i ] = node_handles[i+1];
    }

    return MyMesh( node_handles, coords, quad_handles, quad_connectivity );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

/*
 * This test will repartition via RCB the following quad mesh partitioned on 4
 * processes.

 *-------*-------*-------*-------*
 |       |       |       |       |
 |   0   |   1   |   2   |   3   |
 |       |       |       |       |
 *-------*-------*-------*-------*
 |       |       |       |       |
 |   0   |   1   |   2   |   3   |
 |       |       |       |       |
 *-------*-------*-------*-------*
 |       |       |       |       |
 |   0   |   1   |   2   |   3   |
 |       |       |       |       |
 *-------*-------*-------*-------*
 |       |       |       |       |
 |   0   |   1   |   2   |   3   |
 |       |       |       |       |
 *-------*-------*-------*-------*

 * the resulting repartitioning considering node overlap should be:

 *-------*-------*-------*-------*
 |       |       |       |       |
 |   2   |   2   |   3   |   3   |
 |       |       |       |       |
 *-------*-------*-------*-------*
 |       |       |       |       |
 |   2   |   2   |   3   |   3   |
 |       |       |       |       |
 *-------*-------*-------*-------*
 |       |       |       |       |
 |   0   |   0   |   1   |   1   |
 |       |       |       |       |
 *-------*-------*-------*-------*
 |       |       |       |       |
 |   0   |   0   |   1   |   1   |
 |       |       |       |       |
 *-------*-------*-------*-------*

 */
TEUCHOS_UNIT_TEST( Rendezvous, rendezvous_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // This is a 4 processor test.
    if ( my_size == 4 )
    {
	// Create a bounding box that covers the entire mesh.
	BoundingBox box( -100, -100, -100, 100, 100, 100 );

	// Create a mesh.
	MyMesh my_mesh = buildMyMesh();

	// Create a rendezvous.
	Rendezvous<MyMesh> rendezvous( getDefaultComm<int>(), box );
	rendezvous.build( my_mesh );
	
	// Check the repartitioning with coordinates. Because of the overlap in
	// the rendezvous algorithm, several partitions will contain the proper
	// element, however this is directly checking the RCB algorithm bounding
	// boxes.
	int local_num_quads = 
	    std::distance( my_mesh.quadsBegin(), my_mesh.quadsEnd() );
	int global_num_quads = local_num_quads*my_size;
	int idx;
	int node_dim = 2;
	Teuchos::ArrayRCP<double> coords( node_dim*global_num_quads );
	for ( int j = 0; j < local_num_quads; ++j )
	{
	    for ( int i = 0; i < local_num_quads; ++i )
	    {
		idx = i + local_num_quads*j;
		coords[ idx ] = i + 0.5;
		coords[ global_num_quads + idx ] = j + 0.5;
	    }
	}

	Teuchos::Array<int> destinations = 
	    rendezvous.getRendezvousProcs( coords );
	TEST_ASSERT( destinations.size() == global_num_quads );
	TEST_ASSERT( destinations[0] == 0 );
	TEST_ASSERT( destinations[1] == 0 );
	TEST_ASSERT( destinations[2] == 1 );
	TEST_ASSERT( destinations[3] == 1 );
	TEST_ASSERT( destinations[4] == 0 );
	TEST_ASSERT( destinations[5] == 0 );
	TEST_ASSERT( destinations[6] == 1 );
	TEST_ASSERT( destinations[7] == 1 );
	TEST_ASSERT( destinations[8] == 2 );
	TEST_ASSERT( destinations[9] == 2 );
	TEST_ASSERT( destinations[10] == 3 );
	TEST_ASSERT( destinations[11] == 3 );
	TEST_ASSERT( destinations[12] == 2 );
	TEST_ASSERT( destinations[13] == 2 );
	TEST_ASSERT( destinations[14] == 3 );
	TEST_ASSERT( destinations[15] == 3 );

	// Check the kD-tree with coordinates.
	if ( my_rank == 0 )
	{
	    Teuchos::ArrayRCP<double> 
		target_coords( node_dim*global_num_quads / 4 );

	    target_coords[0] = coords[ 0 ];
	    target_coords[4] = coords[ 0 + global_num_quads ];

	    target_coords[1] = coords[ 1 ];
	    target_coords[5] = coords[ 1 + global_num_quads ];

	    target_coords[2] = coords[ 4 ];
	    target_coords[6] = coords[ 4 + global_num_quads ];

	    target_coords[3] = coords[ 5 ];
	    target_coords[7] = coords[ 5 + global_num_quads ];
	    
	    Teuchos::Array<int> target_elements = 
		rendezvous.getElements( target_coords );

	    TEST_ASSERT( target_elements[0] == 0 );
	    TEST_ASSERT( target_elements[1] == 4 );
	    TEST_ASSERT( target_elements[2] == 1 );
	    TEST_ASSERT( target_elements[3] == 5 );
	}
	if ( my_rank == 1 )
	{
	    Teuchos::ArrayRCP<double> 
		target_coords( node_dim*global_num_quads / 4 );

	    target_coords[0] = coords[ 2 ];
	    target_coords[4] = coords[ 2 + global_num_quads ];

	    target_coords[1] = coords[ 3 ];
	    target_coords[5] = coords[ 3 + global_num_quads ];

	    target_coords[2] = coords[ 6 ];
	    target_coords[6] = coords[ 6 + global_num_quads ];

	    target_coords[3] = coords[ 7 ];
	    target_coords[7] = coords[ 7 + global_num_quads ];
	    
	    Teuchos::Array<int> target_elements = 
		rendezvous.getElements( target_coords );

	    TEST_ASSERT( target_elements[0] == 2 );
	    TEST_ASSERT( target_elements[1] == 3 );
	    TEST_ASSERT( target_elements[2] == 6 );
	    TEST_ASSERT( target_elements[3] == 7 );
	}
	if ( my_rank == 2 )
	{
	    Teuchos::ArrayRCP<double> 
		target_coords( node_dim*global_num_quads / 4 );

	    target_coords[0] = coords[ 8 ];
	    target_coords[4] = coords[ 8 + global_num_quads ];

	    target_coords[1] = coords[ 9 ];
	    target_coords[5] = coords[ 9 + global_num_quads ];

	    target_coords[2] = coords[ 12 ];
	    target_coords[6] = coords[ 12 + global_num_quads ];

	    target_coords[3] = coords[ 13 ];
	    target_coords[7] = coords[ 13 + global_num_quads ];
	    
	    Teuchos::Array<int> target_elements = 
		rendezvous.getElements( target_coords );

	    TEST_ASSERT( target_elements[0] == 8 );
	    TEST_ASSERT( target_elements[1] == 9 );
	    TEST_ASSERT( target_elements[2] == 12 );
	    TEST_ASSERT( target_elements[3] == 13 );
	}
	if ( my_rank == 3 )
	{
	    Teuchos::ArrayRCP<double> 
		target_coords( node_dim*global_num_quads / 4 );

	    target_coords[0] = coords[ 10 ];
	    target_coords[4] = coords[ 10 + global_num_quads ];

	    target_coords[1] = coords[ 11 ];
	    target_coords[5] = coords[ 11 + global_num_quads ];

	    target_coords[2] = coords[ 14 ];
	    target_coords[6] = coords[ 14 + global_num_quads ];

	    target_coords[3] = coords[ 15 ];
	    target_coords[7] = coords[ 15 + global_num_quads ];
	    
	    Teuchos::Array<int> target_elements = 
		rendezvous.getElements( target_coords );

	    TEST_ASSERT( target_elements[0] == 10 );
	    TEST_ASSERT( target_elements[1] == 11 );
	    TEST_ASSERT( target_elements[2] == 14 );
	    TEST_ASSERT( target_elements[3] == 15 );
	}
    }
}
//---------------------------------------------------------------------------//
// end tstRendezvous.cpp
//---------------------------------------------------------------------------//
