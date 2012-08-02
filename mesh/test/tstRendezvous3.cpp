//---------------------------------------------------------------------------//
/*!
 * \file tstRendezvous3.cpp
 * \author Stuart R. Slattery
 * \brief Rendezvous decomposition tests 3.
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

#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_Rendezvous.hpp>
#include <DTK_MeshContainer.hpp>

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
#include <Teuchos_ScalarTraits.hpp>

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
// Mesh create functions.
//---------------------------------------------------------------------------//
DataTransferKit::MeshContainer<int> 
buildTriMesh( int my_rank, int my_size, int edge_length, int elem_offset )
{
    // Make some nodes.
    int num_nodes = edge_length*edge_length;
    int node_dim = 2;
    Teuchos::ArrayRCP<int> node_handles( num_nodes );
    Teuchos::ArrayRCP<double> coords( node_dim*num_nodes );
    int idx;
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length;
	    node_handles[ idx ] = (int) num_nodes*my_rank + idx + elem_offset;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_nodes + idx ] = j;
	}
    }
    
    // Make the triangles. 
    int num_elements = (edge_length-1)*(edge_length-1)*2;
    Teuchos::ArrayRCP<int> tri_handles( num_elements );
    Teuchos::ArrayRCP<int> tri_connectivity( 3*num_elements );
    int elem_idx, node_idx;
    int v0, v1, v2, v3;
    for ( int j = 0; j < (edge_length-1); ++j )
    {
	for ( int i = 0; i < (edge_length-1); ++i )
	{
	    // Indices.
	    node_idx = i + j*edge_length;
	    v0 = node_idx;
	    v1 = node_idx + 1;
	    v2 = node_idx + 1 + edge_length;
	    v3 = node_idx +     edge_length;

	    // Triangle 1.
	    elem_idx = i + j*(edge_length-1);
	    tri_handles[elem_idx] = elem_idx + elem_offset;
	    tri_connectivity[elem_idx]                = node_handles[v0];
	    tri_connectivity[num_elements+elem_idx]   = node_handles[v1];
	    tri_connectivity[2*num_elements+elem_idx] = node_handles[v2];

	    // Triangle 2.
	    elem_idx = i + j*(edge_length-1) + num_elements/2;
	    tri_handles[elem_idx] = elem_idx + elem_offset;
	    tri_connectivity[elem_idx] 	              = node_handles[v2];
	    tri_connectivity[num_elements+elem_idx]   = node_handles[v3];
	    tri_connectivity[2*num_elements+elem_idx] = node_handles[v0];
	}
    }

    Teuchos::ArrayRCP<std::size_t> permutation_list( 3 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return DataTransferKit::MeshContainer<int>( 2, node_handles, coords, 
						DataTransferKit::DTK_TRIANGLE, 3,
						tri_handles, tri_connectivity,
						permutation_list );
}

//---------------------------------------------------------------------------//
DataTransferKit::MeshContainer<int> buildNullTriMesh()
{
    Teuchos::ArrayRCP<int> node_handles(0);
    Teuchos::ArrayRCP<double> coords(0);
    Teuchos::ArrayRCP<int> tri_handles(0);
    Teuchos::ArrayRCP<int> tri_connectivity(0);
    Teuchos::ArrayRCP<std::size_t> permutation_list(3);
    for ( int i = 0; (int) i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return DataTransferKit::MeshContainer<int>( 2, node_handles, coords, 
						DataTransferKit::DTK_TRIANGLE, 3,
						tri_handles, tri_connectivity,
						permutation_list );
}

//---------------------------------------------------------------------------//
DataTransferKit::MeshContainer<int>  
buildQuadMesh( int my_rank, int my_size, int edge_length, int elem_offset )
{
    // Make some nodes.
    int num_nodes = edge_length*edge_length;
    int node_dim = 2;
    Teuchos::ArrayRCP<int> node_handles( num_nodes );
    Teuchos::ArrayRCP<double> coords( node_dim*num_nodes );
    int idx;
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length;
	    node_handles[ idx ] = (int) num_nodes*my_rank + idx + elem_offset;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_nodes + idx ] = j;
	}
    }
    
    // Make the quadrilaterals. 
    int num_elements = (edge_length-1)*(edge_length-1);
    Teuchos::ArrayRCP<int> quad_handles( num_elements );
    Teuchos::ArrayRCP<int> quad_connectivity( 4*num_elements );
    int elem_idx, node_idx;
    for ( int j = 0; j < (edge_length-1); ++j )
    {
	for ( int i = 0; i < (edge_length-1); ++i )
	{
	    node_idx = i + j*edge_length;
	    elem_idx = i + j*(edge_length-1);

	    quad_handles[elem_idx] = elem_idx + elem_offset;

	    quad_connectivity[elem_idx] 
		= node_handles[node_idx];

	    quad_connectivity[num_elements+elem_idx] 
		= node_handles[node_idx+1];

	    quad_connectivity[2*num_elements+elem_idx] 
		= node_handles[node_idx+edge_length+1];

	    quad_connectivity[3*num_elements+elem_idx] 
		= node_handles[node_idx+edge_length];
	}
    }

    Teuchos::ArrayRCP<std::size_t> permutation_list( 4 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return DataTransferKit::MeshContainer<int>( 2, node_handles, coords, 
						DataTransferKit::DTK_QUADRILATERAL, 4,
						quad_handles, quad_connectivity,
						permutation_list );
}

//---------------------------------------------------------------------------//
DataTransferKit::MeshContainer<int> 
buildNullQuadMesh()
{
    Teuchos::ArrayRCP<int> node_handles(0);
    Teuchos::ArrayRCP<double> coords(0);
    Teuchos::ArrayRCP<int> quad_handles(0);
    Teuchos::ArrayRCP<int> quad_connectivity(0);
    Teuchos::ArrayRCP<std::size_t> permutation_list(4);
    for ( int i = 0; (int) i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return DataTransferKit::MeshContainer<int>( 2, node_handles, coords, 
						DataTransferKit::DTK_QUADRILATERAL, 4,
						quad_handles, quad_connectivity,
						permutation_list );
}

//---------------------------------------------------------------------------//
// Global test parameters.
//---------------------------------------------------------------------------//

// number of random points to be generated.
int num_points = 1000;

//---------------------------------------------------------------------------//
// Unit tests.
//---------------------------------------------------------------------------//
// Not all points will be in the mesh in this test.
TEUCHOS_UNIT_TEST( Rendezvous, rendezvous_test3 )
{
    using namespace DataTransferKit;

    // Typedefs.
    typedef MeshContainer<int> MeshType;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Create a bounding box that covers the entire mesh.
    double min = -Teuchos::ScalarTraits<double>::rmax();
    double max = Teuchos::ScalarTraits<double>::rmax();
    BoundingBox box( -100, -100, min, 100, 100, max );

    // Compute element ordinal offsets so we make unique global ordinals.
    int edge_size = 10;
    int tri_offset_1 = 0;
    int quad_offset_1 = tri_offset_1 + (edge_size+1)*(edge_size+1)*2;
    int tri_offset_2 = quad_offset_1 + (edge_size+1)*(edge_size+1);
    int quad_offset_2 = tri_offset_2 + (edge_size+1)*(edge_size+1)*2;

    // Setup source mesh manager.
    Teuchos::ArrayRCP<MeshContainer<int> > mesh_blocks( 2 );
    if ( my_rank == 0 )
    {
	mesh_blocks[0] = 
	    buildTriMesh( my_rank, my_size, edge_size, tri_offset_1 );
	mesh_blocks[1] = buildNullQuadMesh();
    }
    else if ( my_rank == 1 )
    {
	mesh_blocks[0] = buildNullTriMesh();
	mesh_blocks[1] = 
	    buildQuadMesh( my_rank, my_size, edge_size, quad_offset_1 );
    }
    else if ( my_rank == 2 )
    {
	mesh_blocks[0] = 
	    buildTriMesh( my_rank, my_size, edge_size, tri_offset_2 );
	mesh_blocks[1] = buildNullQuadMesh();
    }
    else 
    {
	mesh_blocks[0] = buildNullTriMesh();
	mesh_blocks[1] = 
	    buildQuadMesh( my_rank, my_size, edge_size, quad_offset_2 );
    }
    comm->barrier();

    // Create a mesh manager.
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 2 ) );

    // Create a rendezvous.
    Rendezvous<MeshType> rendezvous( getDefaultComm<int>(), box );
    rendezvous.build( mesh_manager );

    // Give every process a unique block of random numbers;
    std::srand( my_rank*num_points*mesh_manager->dim() );

    // Create some random points.
    int num_rand = num_points*mesh_manager->dim();
    Teuchos::ArrayRCP<double> points( num_rand );
    for ( int i = 0; i < num_points; ++i )
    {
	points[i] = 
	    (my_rank+1)*(edge_size) * (double) std::rand() / RAND_MAX - 1;
	points[i+num_points] = 
	    (edge_size) * (double) std::rand() / RAND_MAX - 0.5;
    }

    // Get the destination procs for the random points and check that they are
    // valid.
    Teuchos::Array<int> destinations = rendezvous.procsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    Teuchos::Array<int> elements, elem_src_procs;
    rendezvous.elementsContainingPoints( points, elements, elem_src_procs );

    // Verify that the points in the mesh were found.
    int num_found = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size*(edge_size-1) < points[i] ||
	     points[num_points + i] < 0.0 || 
	     (edge_size-1) < points[num_points + i] )
	{
	    TEST_ASSERT( elements[i] == -1 );
	}
	else
	{
	    TEST_ASSERT( elements[i] != -1 );
	    ++num_found;
	}
    }
    TEST_ASSERT( num_found > 0 );
}

//---------------------------------------------------------------------------//
// end tstRendezvous3.cpp
//---------------------------------------------------------------------------//
