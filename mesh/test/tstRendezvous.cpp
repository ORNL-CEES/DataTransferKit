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
#include <DTK_MeshManager.hpp>
#include <DTK_MeshTools.hpp>
#include <DTK_MeshContainer.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
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
// Mesh container creation functions.
//---------------------------------------------------------------------------//
// Line segment mesh.
DataTransferKit::MeshContainer<int> buildLineContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some nodes.
    Teuchos::Array<int> node_handles;
    Teuchos::Array<double> coords;

    int node_dim = 1;
    int num_nodes = 2;

    // handles
    for ( int i = 0; i < num_nodes; ++i )
    {
	node_handles.push_back( num_nodes*my_rank + i );
    }

    // x
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 

    // Make the line.
    Teuchos::Array<int> line_handles;
    Teuchos::Array<int> line_connectivity;
    
    // handles
    line_handles.push_back( 12+my_rank );

    // connectivity
    for ( int i = 0; i < num_nodes; ++i )
    {
	line_connectivity.push_back( node_handles[i] );
    }

    Teuchos::ArrayRCP<int> node_handle_array( node_handles.size() );
    std::copy( node_handles.begin(), node_handles.end(), 
	       node_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> line_handle_array( line_handles.size() );
    std::copy( line_handles.begin(), line_handles.end(), 
	       line_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( line_connectivity.size() );
    std::copy( line_connectivity.begin(), line_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<std::size_t> permutation_list( num_nodes );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return MeshContainer<int>( node_dim, node_handle_array, coords_array,
			       DTK_LINE_SEGMENT, num_nodes,
			       line_handle_array, connectivity_array,
			       permutation_list );
}

//---------------------------------------------------------------------------//
// Tri mesh.
DataTransferKit::MeshContainer<int> buildTriContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some nodes.
    Teuchos::Array<int> node_handles;
    Teuchos::Array<double> coords;

    int node_dim = 2;
    int num_nodes = 3;

    // handles
    for ( int i = 0; i < num_nodes; ++i )
    {
	node_handles.push_back( num_nodes*my_rank + i );
    }

    // x
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank+1 ); 

    // y
    coords.push_back( my_rank ); 
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 

    // Make the triahedron.
    Teuchos::Array<int> tri_handles;
    Teuchos::Array<int> tri_connectivity;
    
    // handles
    tri_handles.push_back( 12+my_rank );

    // connectivity
    for ( int i = 0; i < num_nodes; ++i )
    {
	tri_connectivity.push_back( node_handles[i] );
    }
    
    Teuchos::ArrayRCP<int> node_handle_array( node_handles.size() );
    std::copy( node_handles.begin(), node_handles.end(), 
	       node_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> tri_handle_array( tri_handles.size() );
    std::copy( tri_handles.begin(), tri_handles.end(), 
	       tri_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( tri_connectivity.size() );
    std::copy( tri_connectivity.begin(), tri_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<std::size_t> permutation_list( num_nodes );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return MeshContainer<int>( node_dim, node_handle_array, coords_array,
			       DTK_TRIANGLE, num_nodes,
			       tri_handle_array, connectivity_array,
			       permutation_list );
}

//---------------------------------------------------------------------------//
// Quad mesh.
DataTransferKit::MeshContainer<int> buildQuadContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some nodes.
    Teuchos::Array<int> node_handles;
    Teuchos::Array<double> coords;

    int node_dim = 2;
    int num_nodes = 4;

    // handles
    for ( int i = 0; i < num_nodes; ++i )
    {
	node_handles.push_back( num_nodes*my_rank + i );
    }

    // x
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank );

    // y
    coords.push_back( my_rank ); 
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank+1 ); 

    // Make the quadrilateral.
    Teuchos::Array<int> quad_handles;
    Teuchos::Array<int> quad_connectivity;
    
    // handles
    quad_handles.push_back( 12+my_rank );

    // connectivity
    for ( int i = 0; i < num_nodes; ++i )
    {
	quad_connectivity.push_back( node_handles[i] );
    }
    
    Teuchos::ArrayRCP<int> node_handle_array( node_handles.size() );
    std::copy( node_handles.begin(), node_handles.end(), 
	       node_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> quad_handle_array( quad_handles.size() );
    std::copy( quad_handles.begin(), quad_handles.end(), 
	       quad_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( quad_connectivity.size() );
    std::copy( quad_connectivity.begin(), quad_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<std::size_t> permutation_list( num_nodes );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return MeshContainer<int>( node_dim, node_handle_array, coords_array,
			       DTK_QUADRILATERAL, num_nodes,
			       quad_handle_array, connectivity_array,
			       permutation_list );
}

//---------------------------------------------------------------------------//
// Shifted quad mesh.
DataTransferKit::MeshContainer<int> buildShiftedQuadContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some nodes.
    Teuchos::Array<int> node_handles;
    Teuchos::Array<double> coords;

    int node_dim = 2;
    int num_nodes = 4;

    // handles
    for ( int i = 0; i < num_nodes; ++i )
    {
	node_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 );

    // y
    coords.push_back( -1.0 ); 
    coords.push_back( -1.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 

    // Make the quadahedron.
    Teuchos::Array<int> quad_handles;
    Teuchos::Array<int> quad_connectivity;
    
    // handles
    quad_handles.push_back( 9 );

    // connectivity
    for ( int i = 0; i < num_nodes; ++i )
    {
	quad_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<int> node_handle_array( node_handles.size() );
    std::copy( node_handles.begin(), node_handles.end(), 
	       node_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> quad_handle_array( quad_handles.size() );
    std::copy( quad_handles.begin(), quad_handles.end(), 
	       quad_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( quad_connectivity.size() );
    std::copy( quad_connectivity.begin(), quad_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<std::size_t> permutation_list( num_nodes );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return MeshContainer<int>( node_dim, node_handle_array, coords_array,
			       DTK_QUADRILATERAL, num_nodes,
			       quad_handle_array, connectivity_array,
			       permutation_list );
}

//---------------------------------------------------------------------------//
// Tet mesh.
DataTransferKit::MeshContainer<int> buildTetContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some nodes.
    Teuchos::Array<int> node_handles;
    Teuchos::Array<double> coords;

    int node_dim = 3;
    int num_nodes = 4;

    // handles
    for ( int i = 0; i < num_nodes; ++i )
    {
	node_handles.push_back( num_nodes*my_rank + i );
    }

    // x
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank ); 
    coords.push_back( my_rank );

    // y
    coords.push_back( my_rank ); 
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank ); 

    // z
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank+1 );

    // Make the tetahedron.
    Teuchos::Array<int> tet_handles;
    Teuchos::Array<int> tet_connectivity;
    
    // handles
    tet_handles.push_back( 12+my_rank );

    // connectivity
    for ( int i = 0; i < num_nodes; ++i )
    {
	tet_connectivity.push_back( node_handles[i] );
    }
    
    Teuchos::ArrayRCP<int> node_handle_array( node_handles.size() );
    std::copy( node_handles.begin(), node_handles.end(), 
	       node_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> tet_handle_array( tet_handles.size() );
    std::copy( tet_handles.begin(), tet_handles.end(), 
	       tet_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( tet_connectivity.size() );
    std::copy( tet_connectivity.begin(), tet_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<std::size_t> permutation_list( num_nodes );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return MeshContainer<int>( node_dim, node_handle_array, coords_array,
			       DTK_TETRAHEDRON, num_nodes,
			       tet_handle_array, connectivity_array,
			       permutation_list );
}

//---------------------------------------------------------------------------//
// Hex mesh.
DataTransferKit::MeshContainer<int> buildHexContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some nodes.
    Teuchos::Array<int> node_handles;
    Teuchos::Array<double> coords;

    int node_dim = 3;
    int num_nodes = 8;

    // handles
    for ( int i = 0; i < num_nodes; ++i )
    {
	node_handles.push_back( num_nodes*my_rank + i );
    }

    // x
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank ); 

    // y
    coords.push_back( my_rank ); 
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank ); 
    coords.push_back( my_rank );
    coords.push_back( my_rank+1 );
    coords.push_back( my_rank+1 );

    // z
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank+1 );
    coords.push_back( my_rank+1 );
    coords.push_back( my_rank+1 );
    coords.push_back( my_rank+1 );

    // Make the hexahedron.
    Teuchos::Array<int> hex_handles;
    Teuchos::Array<int> hex_connectivity;
    
    // handles
    hex_handles.push_back( 12+my_rank );

    // connectivity
    for ( int i = 0; i < num_nodes; ++i )
    {
	hex_connectivity.push_back( node_handles[i] );
    }
    
    Teuchos::ArrayRCP<int> node_handle_array( node_handles.size() );
    std::copy( node_handles.begin(), node_handles.end(), 
	       node_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> hex_handle_array( hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(), 
	       hex_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<std::size_t> permutation_list( num_nodes );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return MeshContainer<int>( node_dim, node_handle_array, coords_array,
			       DTK_HEXAHEDRON, num_nodes,
			       hex_handle_array, connectivity_array,
			       permutation_list );
}

//---------------------------------------------------------------------------//
// Shifted hex mesh.
DataTransferKit::MeshContainer<int> buildShiftedHexContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some nodes.
    Teuchos::Array<int> node_handles;
    Teuchos::Array<double> coords;

    int node_dim = 3;
    int num_nodes = 8;

    // handles
    for ( int i = 0; i < num_nodes; ++i )
    {
	node_handles.push_back( i );
    }

    // x
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

    // z
    coords.push_back( -1.0 );
    coords.push_back( -1.0 );
    coords.push_back( -1.0 );
    coords.push_back( -1.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );

    // Make the hexahedron.
    Teuchos::Array<int> hex_handles;
    Teuchos::Array<int> hex_connectivity;
    
    // handles
    hex_handles.push_back( 6 );

    // connectivity
    for ( int i = 0; i < num_nodes; ++i )
    {
	hex_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<int> node_handle_array( node_handles.size() );
    std::copy( node_handles.begin(), node_handles.end(), 
	       node_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> hex_handle_array( hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(), 
	       hex_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<std::size_t> permutation_list( num_nodes );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return MeshContainer<int>( node_dim, node_handle_array, coords_array,
			       DTK_HEXAHEDRON, num_nodes,
			       hex_handle_array, connectivity_array,
			       permutation_list );
}

//---------------------------------------------------------------------------//
// Pyramid mesh.
DataTransferKit::MeshContainer<int> buildPyramidContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some nodes.
    Teuchos::Array<int> node_handles;
    Teuchos::Array<double> coords;

    int node_dim = 3;
    int num_nodes = 5;

    // handles
    for ( int i = 0; i < num_nodes; ++i )
    {
	node_handles.push_back( num_nodes*my_rank + i );
    }

    // x
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank );
    coords.push_back( my_rank );

    // y
    coords.push_back( my_rank ); 
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank ); 

    // z
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank+1 );

    // Make the pyramidahedron.
    Teuchos::Array<int> pyramid_handles;
    Teuchos::Array<int> pyramid_connectivity;
    
    // handles
    pyramid_handles.push_back( 12+my_rank );

    // connectivity
    for ( int i = 0; i < num_nodes; ++i )
    {
	pyramid_connectivity.push_back( node_handles[i] );
    }
    
    Teuchos::ArrayRCP<int> node_handle_array( node_handles.size() );
    std::copy( node_handles.begin(), node_handles.end(), 
	       node_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> pyramid_handle_array( pyramid_handles.size() );
    std::copy( pyramid_handles.begin(), pyramid_handles.end(), 
	       pyramid_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( pyramid_connectivity.size() );
    std::copy( pyramid_connectivity.begin(), pyramid_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<std::size_t> permutation_list( num_nodes );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return MeshContainer<int>( node_dim, node_handle_array, coords_array,
			       DTK_PYRAMID, num_nodes,
			       pyramid_handle_array, connectivity_array,
			       permutation_list );
}

//---------------------------------------------------------------------------//
// Shifted pyramid mesh.
DataTransferKit::MeshContainer<int> buildShiftedPyramidContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some nodes.
    Teuchos::Array<int> node_handles;
    Teuchos::Array<double> coords;

    int node_dim = 3;
    int num_nodes = 5;

    // handles
    for ( int i = 0; i < num_nodes; ++i )
    {
	node_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 );
    coords.push_back( 0.0 );

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 

    // z
    coords.push_back( -1.0 );
    coords.push_back( -1.0 );
    coords.push_back( -1.0 );
    coords.push_back( -1.0 );
    coords.push_back( -2.0 );

    // Make the pyramidahedron.
    Teuchos::Array<int> pyramid_handles;
    Teuchos::Array<int> pyramid_connectivity;
    
    // handles
    pyramid_handles.push_back( 89 );

    // connectivity
    for ( int i = 0; i < num_nodes; ++i )
    {
	pyramid_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<int> node_handle_array( node_handles.size() );
    std::copy( node_handles.begin(), node_handles.end(), 
	       node_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> pyramid_handle_array( pyramid_handles.size() );
    std::copy( pyramid_handles.begin(), pyramid_handles.end(), 
	       pyramid_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( pyramid_connectivity.size() );
    std::copy( pyramid_connectivity.begin(), pyramid_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<std::size_t> permutation_list( num_nodes );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return MeshContainer<int>( node_dim, node_handle_array, coords_array,
			       DTK_PYRAMID, num_nodes,
			       pyramid_handle_array, connectivity_array,
			       permutation_list );
}

//---------------------------------------------------------------------------//
// Global Test Variables.
//---------------------------------------------------------------------------//

int num_points = 1000;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Line mesh.
TEUCHOS_UNIT_TEST( MeshContainer, line_rendezvous_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Create a bounding box that covers the entire mesh.
    double min = -Teuchos::ScalarTraits<double>::rmax();
    double max = Teuchos::ScalarTraits<double>::rmax();
    BoundingBox box( -100, min, min, 100, max, max );

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP< MeshType > mesh_blocks( 1 );
    mesh_blocks[0] = buildLineContainer( my_rank );

    // Create a mesh manager.
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 1 ) );

    // Create a rendezvous.
    Rendezvous<MeshType> rendezvous( getDefaultComm<int>(), box );
    rendezvous.build( mesh_manager );

    // Give every process a unique block of random numbers;
    std::srand( my_rank*num_points*mesh_manager->dim() );

    // Create some random points.
    Teuchos::ArrayRCP<double> points(num_points);
    for ( int i = 0; i < num_points; ++i )
    {
	points[i] = (my_size+1) * (double) std::rand() / RAND_MAX - 0.5;
    }

    // Get the destination procs for the random points and check that its in a
    // valid range.
    Teuchos::Array<int> destinations = rendezvous.getRendezvousProcs( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points and check
    // that they are found in the correct element.
    Teuchos::Array<int> elements = rendezvous.elementsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] )
	{
	    TEST_ASSERT( elements[i] == -1 );
	}
	else if ( elements[i] != -1 )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Tri mesh.
TEUCHOS_UNIT_TEST( MeshContainer, tri_rendezvous_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Create a bounding box that covers the entire mesh.
    BoundingBox box( -100, -100, -100, 100, 100, 100 );

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP< MeshType > mesh_blocks( 1 );
    mesh_blocks[0] = buildTriContainer( my_rank );

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
    for ( int i = 0; i < num_rand; ++i )
    {
	points[i] = (my_size+1) * (double) std::rand() / RAND_MAX - 0.5;
    }

    // Get the destination procs for the random points.
    Teuchos::Array<int> destinations = rendezvous.getRendezvousProcs( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    Teuchos::Array<int> elements = rendezvous.elementsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] ||
	     points[num_points + i] < 0.0 || my_size < points[num_points + i] )
	{
	    TEST_ASSERT( elements[i] == -1 );
	}
	else if ( elements[i] != -1 )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Quad mesh.
TEUCHOS_UNIT_TEST( MeshContainer, quad_rendezvous_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Create a bounding box that covers the entire mesh.
    BoundingBox box( -100, -100, -100, 100, 100, 100 );

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP< MeshType > mesh_blocks( 1 );
    mesh_blocks[0] = buildQuadContainer( my_rank );

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
    for ( int i = 0; i < num_rand; ++i )
    {
	points[i] = (my_size+1) * (double) std::rand() / RAND_MAX - 0.5;
    }

    // Get the destination procs for the random points.
    Teuchos::Array<int> destinations = rendezvous.getRendezvousProcs( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    Teuchos::Array<int> elements = rendezvous.elementsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] ||
	     points[num_points + i] < 0.0 || my_size < points[num_points + i] )
	{
	    TEST_ASSERT( elements[i] == -1 );
	}
	else if ( elements[i] != -1 )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Tet mesh.
TEUCHOS_UNIT_TEST( MeshContainer, tet_rendezvous_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Create a bounding box that covers the entire mesh.
    BoundingBox box( -100, -100, -100, 100, 100, 100 );

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP< MeshType > mesh_blocks( 1 );
    mesh_blocks[0] = buildTetContainer( my_rank );

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
    for ( int i = 0; i < num_rand; ++i )
    {
	points[i] = (my_size+1) * (double) std::rand() / RAND_MAX - 0.5;
    }

    // Get the destination procs for the random points.
    Teuchos::Array<int> destinations = rendezvous.getRendezvousProcs( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    Teuchos::Array<int> elements = rendezvous.elementsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] ||
	     points[num_points + i] < 0.0 || my_size < points[num_points + i] ||
	     points[2*num_points + i] < 0.0 || my_size < points[2*num_points + i] )
	{
	    TEST_ASSERT( elements[i] == -1 );
	}
	else if ( elements[i] != -1 )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Hex mesh.
TEUCHOS_UNIT_TEST( MeshContainer, hex_rendezvous_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Create a bounding box that covers the entire mesh.
    BoundingBox box( -100, -100, -100, 100, 100, 100 );

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP< MeshType > mesh_blocks( 1 );
    mesh_blocks[0] = buildHexContainer( my_rank );

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
    for ( int i = 0; i < num_rand; ++i )
    {
	points[i] = (my_size+1) * (double) std::rand() / RAND_MAX - 0.5;
    }

    // Get the destination procs for the random points.
    Teuchos::Array<int> destinations = rendezvous.getRendezvousProcs( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    Teuchos::Array<int> elements = rendezvous.elementsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] ||
	     points[num_points + i] < 0.0 || my_size < points[num_points + i] ||
	     points[2*num_points + i] < 0.0 || my_size < points[2*num_points + i] )
	{
	    TEST_ASSERT( elements[i] == -1 );
	}
	else if ( elements[i] != -1 )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Pyramid mesh.
TEUCHOS_UNIT_TEST( MeshContainer, pyramid_rendezvous_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Create a bounding box that covers the entire mesh.
    BoundingBox box( -100, -100, -100, 100, 100, 100 );

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP< MeshType > mesh_blocks( 1 );
    mesh_blocks[0] = buildPyramidContainer( my_rank );

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
    for ( int i = 0; i < num_rand; ++i )
    {
	points[i] = (my_size+1) * (double) std::rand() / RAND_MAX - 0.5;
    }

    // Get the destination procs for the random points.
    Teuchos::Array<int> destinations = rendezvous.getRendezvousProcs( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    Teuchos::Array<int> elements = rendezvous.elementsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] ||
	     points[num_points + i] < 0.0 || my_size < points[num_points + i] ||
	     points[2*num_points + i] < 0.0 || my_size < points[2*num_points + i] )
	{
	    TEST_ASSERT( elements[i] == -1 );
	}
	else if ( elements[i] != -1 )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) );
	}
    }
}

// //---------------------------------------------------------------------------//
// // 2d hybrid test.
// TEUCHOS_UNIT_TEST( MeshContainer, 2d_hybrid_rendezvous_test )
// {
//     using namespace DataTransferKit;

//     // Create a mesh container.
//     typedef MeshContainer<int> MeshType;
//     typedef MeshTraits< MeshType > MT;
//     typedef MeshTools< MeshType > Tools;
//     Teuchos::ArrayRCP< MeshType > mesh_blocks( 2 );
//     mesh_blocks[0] = buildTriContainer( my_rank );
//     mesh_blocks[1] = buildShiftedQuadContainer( my_rank );

//     // Create a mesh manager.
//     MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 2 );

//     // Create a rendezvous.
//     Rendezvous<MeshType> rendezvous( getDefaultComm<int>(), box );
//     rendezvous.build( mesh_manager );

    // Give every process a unique block of random numbers;
//    std::srand( my_rank*num_points*mesh_manager->dim() );


//     // Search the rendezvous decomposition for some random points.
//     Teuchos::Array<double> point(2);
//     int ordinal = 0;
//     double tol = 1.0e-8;
//     for ( int i = 0; i < num_points; ++i )
//     {
// 	ordinal = 0;
// 	point[0] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
// 	point[1] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;

// 	// We can end up either in the quad or tri on their boundary.
// 	if ( 0.0 <= point[0] && point[0] <= 1.0 &&
// 	     -tol <= point[1] && point[1] <= tol )
// 	{
// 	    TEST_ASSERT( kd_tree.findPoint( point, ordinal ) );
// 	    TEST_ASSERT( ordinal == 12 || ordinal == 9 );
// 	}
// 	// In the tri.
// 	else if ( 0.0 <= point[0] && point[0] <= 1.0 &&
// 		  0.0 < point[1] && point[1] <= point[0] )
// 	{
// 	    TEST_ASSERT( kd_tree.findPoint( point, ordinal ) );
// 	    TEST_ASSERT( ordinal == 12 );
// 	}
// 	// In the quad.
// 	else if ( 0.0 <= point[0] && point[0] <= 1.0 &&
// 		  -1.0 <= point[1] && point[1] < 0.0 )
// 	{
// 	    TEST_ASSERT( kd_tree.findPoint( point, ordinal ) );
// 	    TEST_ASSERT( ordinal == 9 );
// 	}
// 	// Neither
// 	else
// 	{
// 	    TEST_ASSERT( !kd_tree.findPoint( point, ordinal ) );
// 	}
//     }
// }

// //---------------------------------------------------------------------------//
// // 3d hybrid test.
// TEUCHOS_UNIT_TEST( MeshContainer, 3d_hybrid_rendezvous_test )
// {
//     using namespace DataTransferKit;

//     int my_rank = getDefaultComm<int>()->getRank();
//     int my_size = getDefaultComm<int>()->getSize();

//     // Create a bounding box that covers the entire mesh.
//     BoundingBox box( -100, -100, -100, 100, 100, 100 );

//     // Create a mesh container.
//     typedef MeshContainer<int> MeshType;
//     typedef MeshTraits< MeshType > MT;
//     typedef MeshTools< MeshType > Tools;
//     Teuchos::ArrayRCP< MeshType > mesh_blocks( 3 );
//     mesh_blocks[0] = buildTetContainer( my_rank );
//     mesh_blocks[1] = buildShiftedPyramidContainer( my_rank );
//     mesh_blocks[2] = buildShiftedHexContainer( my_rank );

//     // Create a mesh manager.
//     MeshManager<MeshType> mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );

//     // Create a rendezvous.
//     Rendezvous<MeshType> rendezvous( getDefaultComm<int>(), box );
//     rendezvous.build( mesh_manager );

    // Give every process a unique block of random numbers;
// std::srand( my_rank*num_points*mesh_manager->dim() );


//     // Search the rendezvous decomposition for some random points.
//     double tol = 1.0e-8;
//     Teuchos::Array<double> point(3);
//     int ordinal = 0;
//     for ( int i = 0; i < num_points; ++i )
//     {
// 	ordinal = 0;
// 	point[0] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
// 	point[1] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
// 	point[2] = 4.0 * (double) std::rand() / RAND_MAX - 2.5;
    	
// 	// Hex/Pyramid boundary.
// 	if ( 0.0 <= point[0] && point[0] <= 1.0 &&
// 	     0.0 <= point[1] && point[1] <= 1.0 &&
// 	     -1.0-tol <= point[2] && point[2] <= -1.0+tol )
// 	{
// 	    TEST_ASSERT( kd_tree.findPoint( point, ordinal ) );
// 	    TEST_ASSERT( ordinal == 6 || ordinal == 89 );
// 	}

// 	// Hex/Tet boundary.
// 	else if ( 0.0 <= point[0] && point[0] <= 1.0 &&
// 		  0.0 <= point[1] && point[1] <= 1.0 &&
// 		  -tol <= point[2] && point[2] <= tol )
// 	{
// 	    TEST_ASSERT( kd_tree.findPoint( point, ordinal ) );
// 	    TEST_ASSERT( ordinal == 12 || ordinal == 6 );
// 	}

// 	// Tet
// 	else if ( std::max( std::max(-point[0],-point[1]),
// 			    std::max(-point[2], point[0]+point[1]+point[2]-1) )
// 		  < tol )
// 	{
// 	    TEST_ASSERT( kd_tree.findPoint( point, ordinal ) );
// 	    TEST_ASSERT( ordinal == 12 );
// 	}

// 	// Hex
// 	else if ( 0.0 <= point[0] && point[0] <= 1.0 &&
// 		  0.0 <= point[1] && point[1] <= 1.0 && 
// 		  -1.0 <= point[2] && point[2] <= 0.0 )
// 	{
// 	    TEST_ASSERT( kd_tree.findPoint( point, ordinal ) );
// 	    TEST_ASSERT( ordinal == 6 );
// 	}
// 	// Pyramid
// 	else if( 0.0 <= point[0] && point[0] <= 2.0+point[2] &&
// 		 0.0 <= point[1] && point[1] <= 2.0+point[2] && 
// 		 -2.0 <= point[2] && point[2] <= -1.0 )
// 	{
// 	    TEST_ASSERT( kd_tree.findPoint( point, ordinal ) );
// 	    TEST_ASSERT( ordinal == 89 );
// 	}

// 	// None
// 	else
// 	{
// 	    TEST_ASSERT( !kd_tree.findPoint( point, ordinal ) );
// 	}
//     }
// }

//---------------------------------------------------------------------------//
// end tstRendezvous.cpp
//---------------------------------------------------------------------------//
