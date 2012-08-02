//---------------------------------------------------------------------------//
/*!
 * \file tstRendezvous2.cpp
 * \author Stuart R. Slattery
 * \brief Rendezvous decomposition tests 2.
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
buildTetMesh( int my_rank, int my_size, int edge_length, int elem_offset )
{
    // Make some nodes.
    int num_nodes = edge_length*edge_length*2;
    int node_dim = 3;
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
	    coords[ 2*num_nodes + idx ] = 0.0;
	}
    }
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length + num_nodes / 2;
	    node_handles[ idx ] = (int) num_nodes*my_rank + idx + elem_offset;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_nodes + idx ] = j;
	    coords[ 2*num_nodes + idx ] = 1.0;
	}
    }
    
    // Make the tetrahedrons. 
    int num_elements = (edge_length-1)*(edge_length-1)*5;
    Teuchos::ArrayRCP<int> tet_handles( num_elements );
    Teuchos::ArrayRCP<int> tet_connectivity( 4*num_elements );
    int elem_idx, node_idx;
    int v0, v1, v2, v3, v4, v5, v6, v7;
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
	    v4 = node_idx +                   num_nodes/2;
	    v5 = node_idx + 1 +               num_nodes/2;
	    v6 = node_idx + 1 + edge_length + num_nodes/2;
	    v7 = node_idx +     edge_length + num_nodes/2; 

	    // Tetrahedron 1.
	    elem_idx = i + j*(edge_length-1);
	    tet_handles[elem_idx] = elem_idx + elem_offset;
	    tet_connectivity[elem_idx]                = node_handles[v0];
	    tet_connectivity[num_elements+elem_idx]   = node_handles[v1];
	    tet_connectivity[2*num_elements+elem_idx] = node_handles[v3];
	    tet_connectivity[3*num_elements+elem_idx] = node_handles[v4];

	    // Tetrahedron 2.
	    elem_idx = i + j*(edge_length-1) + num_elements/5;
	    tet_handles[elem_idx] = elem_idx + elem_offset;
	    tet_connectivity[elem_idx] 	              = node_handles[v1];
	    tet_connectivity[num_elements+elem_idx]   = node_handles[v2];
	    tet_connectivity[2*num_elements+elem_idx] = node_handles[v3];
	    tet_connectivity[3*num_elements+elem_idx] = node_handles[v6];

	    // Tetrahedron 3.
	    elem_idx = i + j*(edge_length-1) + 2*num_elements/5;
	    tet_handles[elem_idx] = elem_idx + elem_offset;
	    tet_connectivity[elem_idx] 	              = node_handles[v6];
	    tet_connectivity[num_elements+elem_idx]   = node_handles[v5];
	    tet_connectivity[2*num_elements+elem_idx] = node_handles[v4];
	    tet_connectivity[3*num_elements+elem_idx] = node_handles[v1];

	    // Tetrahedron 4.
	    elem_idx = i + j*(edge_length-1) + 3*num_elements/5;
	    tet_handles[elem_idx] = elem_idx + elem_offset;
	    tet_connectivity[elem_idx]   	      = node_handles[v4];
	    tet_connectivity[num_elements+elem_idx]   = node_handles[v7];
	    tet_connectivity[2*num_elements+elem_idx] = node_handles[v6];
	    tet_connectivity[3*num_elements+elem_idx] = node_handles[v3];

	    // Tetrahedron 5.
	    elem_idx = i + j*(edge_length-1) + 4*num_elements/5;
	    tet_handles[elem_idx] = elem_idx + elem_offset;
	    tet_connectivity[elem_idx] 	              = node_handles[v3];
	    tet_connectivity[num_elements+elem_idx]   = node_handles[v1];
	    tet_connectivity[2*num_elements+elem_idx] = node_handles[v6];
	    tet_connectivity[3*num_elements+elem_idx] = node_handles[v4];
	}
    }

    Teuchos::ArrayRCP<std::size_t> permutation_list( 4 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return DataTransferKit::MeshContainer<int>( 3, node_handles, coords, 
						DataTransferKit::DTK_TETRAHEDRON, 4,
						tet_handles, tet_connectivity,
						permutation_list );
}

//---------------------------------------------------------------------------//
DataTransferKit::MeshContainer<int> buildNullTetMesh()
{
    Teuchos::ArrayRCP<int> node_handles(0);
    Teuchos::ArrayRCP<double> coords(0);
    Teuchos::ArrayRCP<int> tet_handles(0);
    Teuchos::ArrayRCP<int> tet_connectivity(0);
    Teuchos::ArrayRCP<std::size_t> permutation_list(4);
    for ( int i = 0; (int) i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return DataTransferKit::MeshContainer<int>( 3, node_handles, coords, 
						DataTransferKit::DTK_TETRAHEDRON, 4,
						tet_handles, tet_connectivity,
						permutation_list );
}

//---------------------------------------------------------------------------//
DataTransferKit::MeshContainer<int>  
buildHexMesh( int my_rank, int my_size, int edge_length, int elem_offset )
{
    // Make some nodes.
    int num_nodes = edge_length*edge_length*2;
    int node_dim = 3;
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
	    coords[ 2*num_nodes + idx ] = 0.0;
	}
    }
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length + num_nodes / 2;
	    node_handles[ idx ] = (int) num_nodes*my_rank + idx + elem_offset;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_nodes + idx ] = j;
	    coords[ 2*num_nodes + idx ] = 1.0;
	}
    }
    
    // Make the hexahedrons. 
    int num_elements = (edge_length-1)*(edge_length-1);
    Teuchos::ArrayRCP<int> hex_handles( num_elements );
    Teuchos::ArrayRCP<int> hex_connectivity( 8*num_elements );
    int elem_idx, node_idx;
    for ( int j = 0; j < (edge_length-1); ++j )
    {
	for ( int i = 0; i < (edge_length-1); ++i )
	{
	    node_idx = i + j*edge_length;
	    elem_idx = i + j*(edge_length-1);

	    hex_handles[elem_idx] = elem_idx + elem_offset;

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

    Teuchos::ArrayRCP<std::size_t> permutation_list( 8 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return DataTransferKit::MeshContainer<int>( 3, node_handles, coords, 
						DataTransferKit::DTK_HEXAHEDRON, 8,
						hex_handles, hex_connectivity,
						permutation_list );
}

//---------------------------------------------------------------------------//
DataTransferKit::MeshContainer<int> 
buildNullHexMesh()
{
    Teuchos::ArrayRCP<int> node_handles(0);
    Teuchos::ArrayRCP<double> coords(0);
    Teuchos::ArrayRCP<int> hex_handles(0);
    Teuchos::ArrayRCP<int> hex_connectivity(0);
    Teuchos::ArrayRCP<std::size_t> permutation_list(8);
    for ( int i = 0; (int) i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return DataTransferKit::MeshContainer<int>( 3, node_handles, coords, 
						DataTransferKit::DTK_HEXAHEDRON, 8,
						hex_handles, hex_connectivity,
						permutation_list );
}

//---------------------------------------------------------------------------//
DataTransferKit::MeshContainer<int>  
buildPyramidMesh( int my_rank, int my_size, int edge_length, int elem_offset )
{
    // Make some nodes.
    int num_nodes = edge_length*edge_length*2 + (edge_length-1)*(edge_length-1);
    int node_dim = 3;
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
	    coords[ 2*num_nodes + idx ] = 0.0;
	}
    }
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length + edge_length*edge_length;
	    node_handles[ idx ] = (int) num_nodes*my_rank + idx + elem_offset;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_nodes + idx ] = j;
	    coords[ 2*num_nodes + idx ] = 1.0;
	}
    }
    for ( int j = 0; j < edge_length-1; ++j )
    {
	for ( int i = 0; i < edge_length-1; ++i )
	{
	    idx = i + j*(edge_length-1) + edge_length*edge_length*2;
	    node_handles[ idx ] = (int) num_nodes*my_rank + idx + elem_offset;
	    coords[ idx ] = i + my_rank*(edge_length-1) + 0.5;
	    coords[ num_nodes + idx ] = j + 0.5;
	    coords[ 2*num_nodes + idx ] = 0.5;
	}
    }
    
    // Make the pyramids. 
    int num_elements = (edge_length-1)*(edge_length-1)*6;
    Teuchos::ArrayRCP<int> pyr_handles( num_elements );
    Teuchos::ArrayRCP<int> pyr_connectivity( 5*num_elements );
    int elem_idx, node_idx;
    int v0, v1, v2, v3, v4, v5, v6, v7, v8;
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
	    v4 = node_idx +                   edge_length*edge_length;
	    v5 = node_idx + 1 +               edge_length*edge_length;
	    v6 = node_idx + 1 + edge_length + edge_length*edge_length;
	    v7 = node_idx +     edge_length + edge_length*edge_length;
	    v8 = i + j*(edge_length-1) + edge_length*edge_length*2;

	    // Pyramid 1.
	    elem_idx = i + j*(edge_length-1);
	    pyr_handles[elem_idx] = elem_idx + elem_offset;
	    pyr_connectivity[elem_idx]                = node_handles[v0];
	    pyr_connectivity[num_elements+elem_idx]   = node_handles[v1];
	    pyr_connectivity[2*num_elements+elem_idx] = node_handles[v2];
	    pyr_connectivity[3*num_elements+elem_idx] = node_handles[v3];
	    pyr_connectivity[4*num_elements+elem_idx] = node_handles[v8];

	    // Pyramid 2.
	    elem_idx = i + j*(edge_length-1) + num_elements/6;
	    pyr_handles[elem_idx] = elem_idx + elem_offset;
	    pyr_connectivity[elem_idx] 	              = node_handles[v1];
	    pyr_connectivity[num_elements+elem_idx]   = node_handles[v5];
	    pyr_connectivity[2*num_elements+elem_idx] = node_handles[v6];
	    pyr_connectivity[3*num_elements+elem_idx] = node_handles[v2];
	    pyr_connectivity[4*num_elements+elem_idx] = node_handles[v8];

	    // Pyramid 3.
	    elem_idx = i + j*(edge_length-1) + 2*num_elements/6;
	    pyr_handles[elem_idx] = elem_idx + elem_offset;
	    pyr_connectivity[elem_idx] 	              = node_handles[v2];
	    pyr_connectivity[num_elements+elem_idx]   = node_handles[v6];
	    pyr_connectivity[2*num_elements+elem_idx] = node_handles[v7];
	    pyr_connectivity[3*num_elements+elem_idx] = node_handles[v3];
	    pyr_connectivity[4*num_elements+elem_idx] = node_handles[v8];

	    // Pyramid 4.
	    elem_idx = i + j*(edge_length-1) + 3*num_elements/6;
	    pyr_handles[elem_idx] = elem_idx + elem_offset;
	    pyr_connectivity[elem_idx]   	      = node_handles[v4];
	    pyr_connectivity[num_elements+elem_idx]   = node_handles[v0];
	    pyr_connectivity[2*num_elements+elem_idx] = node_handles[v3];
	    pyr_connectivity[3*num_elements+elem_idx] = node_handles[v7];
	    pyr_connectivity[4*num_elements+elem_idx] = node_handles[v8];

	    // Pyramid 5.
	    elem_idx = i + j*(edge_length-1) + 4*num_elements/6;
	    pyr_handles[elem_idx] = elem_idx + elem_offset;
	    pyr_connectivity[elem_idx]   	      = node_handles[v4];
	    pyr_connectivity[num_elements+elem_idx]   = node_handles[v5];
	    pyr_connectivity[2*num_elements+elem_idx] = node_handles[v1];
	    pyr_connectivity[3*num_elements+elem_idx] = node_handles[v0];
	    pyr_connectivity[4*num_elements+elem_idx] = node_handles[v8];

	    // Pyramid 6.
	    elem_idx = i + j*(edge_length-1) + 5*num_elements/6;
	    pyr_handles[elem_idx] = elem_idx + elem_offset;
	    pyr_connectivity[elem_idx]   	      = node_handles[v4];
	    pyr_connectivity[num_elements+elem_idx]   = node_handles[v7];
	    pyr_connectivity[2*num_elements+elem_idx] = node_handles[v6];
	    pyr_connectivity[3*num_elements+elem_idx] = node_handles[v5];
	    pyr_connectivity[4*num_elements+elem_idx] = node_handles[v8];
	}
    }

    Teuchos::ArrayRCP<std::size_t> permutation_list( 5 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return DataTransferKit::MeshContainer<int>( 3, node_handles, coords, 
						DataTransferKit::DTK_PYRAMID, 5,
						pyr_handles, pyr_connectivity,
						permutation_list );
}

//---------------------------------------------------------------------------//
DataTransferKit::MeshContainer<int> buildNullPyramidMesh()
{
    Teuchos::ArrayRCP<int> node_handles(0);
    Teuchos::ArrayRCP<double> coords(0);
    Teuchos::ArrayRCP<int> pyramid_handles(0);
    Teuchos::ArrayRCP<int> pyramid_connectivity(0);
    Teuchos::ArrayRCP<std::size_t> permutation_list(5);
    for ( int i = 0; (int) i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return DataTransferKit::MeshContainer<int>( 3, node_handles, coords, 
						DataTransferKit::DTK_PYRAMID, 5,
						pyramid_handles, pyramid_connectivity,
						permutation_list );
}

//---------------------------------------------------------------------------//
DataTransferKit::MeshContainer<int>  
buildWedgeMesh( int my_rank, int my_size, int edge_length, int elem_offset )
{
    // Make some nodes.
    int num_nodes = edge_length*edge_length*2;
    int node_dim = 3;
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
	    coords[ 2*num_nodes + idx ] = 0.0;
	}
    }
    for ( int j = 0; j < edge_length; ++j )
    {
	for ( int i = 0; i < edge_length; ++i )
	{
	    idx = i + j*edge_length + edge_length*edge_length;
	    node_handles[ idx ] = (int) num_nodes*my_rank + idx + elem_offset;
	    coords[ idx ] = i + my_rank*(edge_length-1);
	    coords[ num_nodes + idx ] = j;
	    coords[ 2*num_nodes + idx ] = 1.0;
	}
    }
    
    // Make the wedges. 
    int num_elements = (edge_length-1)*(edge_length-1)*2;
    Teuchos::ArrayRCP<int> wedge_handles( num_elements );
    Teuchos::ArrayRCP<int> wedge_connectivity( 6*num_elements );
    int elem_idx, node_idx;
    int v0, v1, v2, v3, v4, v5, v6, v7;
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
	    v4 = node_idx +                   edge_length*edge_length;
	    v5 = node_idx + 1 +               edge_length*edge_length;
	    v6 = node_idx + 1 + edge_length + edge_length*edge_length;
	    v7 = node_idx +     edge_length + edge_length*edge_length;

	    // Wedge 1.
	    elem_idx = i + j*(edge_length-1);
	    wedge_handles[elem_idx] = elem_idx + elem_offset;
	    wedge_connectivity[elem_idx]                = node_handles[v0];
	    wedge_connectivity[num_elements+elem_idx]   = node_handles[v4];
	    wedge_connectivity[2*num_elements+elem_idx] = node_handles[v1];
	    wedge_connectivity[3*num_elements+elem_idx] = node_handles[v3];
	    wedge_connectivity[4*num_elements+elem_idx] = node_handles[v7];
	    wedge_connectivity[5*num_elements+elem_idx] = node_handles[v2];

	    // Wedge 2.
	    elem_idx = i + j*(edge_length-1) + num_elements/2;
	    wedge_handles[elem_idx] = elem_idx + elem_offset;
	    wedge_connectivity[elem_idx] 	        = node_handles[v1];
	    wedge_connectivity[num_elements+elem_idx]   = node_handles[v4];
	    wedge_connectivity[2*num_elements+elem_idx] = node_handles[v5];
	    wedge_connectivity[3*num_elements+elem_idx] = node_handles[v2];
	    wedge_connectivity[4*num_elements+elem_idx] = node_handles[v7];
	    wedge_connectivity[5*num_elements+elem_idx] = node_handles[v6];
	}
    }

    Teuchos::ArrayRCP<std::size_t> permutation_list( 6 );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return DataTransferKit::MeshContainer<int>( 3, node_handles, coords, 
						DataTransferKit::DTK_WEDGE, 6,
						wedge_handles, wedge_connectivity,
						permutation_list );
}

//---------------------------------------------------------------------------//
DataTransferKit::MeshContainer<int> buildNullWedgeMesh()
{
    Teuchos::ArrayRCP<int> node_handles(0);
    Teuchos::ArrayRCP<double> coords(0);
    Teuchos::ArrayRCP<int> wedge_handles(0);
    Teuchos::ArrayRCP<int> wedge_connectivity(0);
    Teuchos::ArrayRCP<std::size_t> permutation_list(6);
    for ( int i = 0; (int) i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }

    return DataTransferKit::MeshContainer<int>( 3, node_handles, coords, 
						DataTransferKit::DTK_WEDGE, 6,
						wedge_handles, wedge_connectivity,
						permutation_list );
}

//---------------------------------------------------------------------------//
// Global test parameters.
//---------------------------------------------------------------------------//

// number of random points to be generated.
int num_points = 500;

//---------------------------------------------------------------------------//
// Unit tests.
//---------------------------------------------------------------------------//
// Not all points will be in the mesh in this test.
TEUCHOS_UNIT_TEST( Rendezvous, rendezvous_test2 )
{
    using namespace DataTransferKit;

    // Typedefs.
    typedef MeshContainer<int> MeshType;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Create a bounding box that covers the entire mesh.
    BoundingBox box( -100, -100, -100, 100, 100, 100 );

    // Compute element ordinal offsets so we make unique global ordinals.
    int edge_size = 10;
    int tet_offset = 0;
    int hex_offset = tet_offset + (edge_size+1)*(edge_size+1)*5;
    int pyramid_offset = hex_offset + (edge_size+1)*(edge_size+1);
    int wedge_offset = pyramid_offset + (edge_size+1)*(edge_size+1)*6;

    // Setup source mesh manager.
    Teuchos::ArrayRCP<MeshContainer<int> > mesh_blocks( 4 );
    if ( my_rank == 0 )
    {
	mesh_blocks[0] = 
	    buildTetMesh( my_rank, my_size, edge_size, tet_offset );
	mesh_blocks[1] = buildNullHexMesh();
	mesh_blocks[2] = buildNullPyramidMesh();
	mesh_blocks[3] = buildNullWedgeMesh();
    }
    else if ( my_rank == 1 )
    {
	mesh_blocks[0] = buildNullTetMesh();
	mesh_blocks[1] = 
	    buildHexMesh( my_rank, my_size, edge_size, hex_offset );
	mesh_blocks[2] = buildNullPyramidMesh();
	mesh_blocks[3] = buildNullWedgeMesh();
    }
    else if ( my_rank == 2 )
    {
	mesh_blocks[0] = buildNullTetMesh();
	mesh_blocks[1] = buildNullHexMesh();
	mesh_blocks[2] = 
	    buildPyramidMesh( my_rank, my_size, edge_size, pyramid_offset );
	mesh_blocks[3] = buildNullWedgeMesh();
    }
    else 
    {
	mesh_blocks[0] = buildNullTetMesh();
	mesh_blocks[1] = buildNullHexMesh();
	mesh_blocks[2] = buildNullPyramidMesh();
	mesh_blocks[3] = 
	    buildWedgeMesh( my_rank, my_size, edge_size, wedge_offset );
    }
    comm->barrier();

    // Create a mesh manager.
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 3 ) );

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
	points[i+2*num_points] = 
	    1.5 * (double) std::rand() / RAND_MAX - 0.25;
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
	     (edge_size-1) < points[num_points + i] ||
	     points[2*num_points + i] < 0.0 || 1.0 < points[2*num_points + i] )
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
// end tstRendezvous2.cpp
//---------------------------------------------------------------------------//
