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
#include <limits>

#include <DTK_Rendezvous.hpp>
#include <DTK_BoundingBox.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshManager.hpp>
#include <DTK_MeshTools.hpp>
#include <DTK_MeshContainer.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
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
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildLineContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 1;
    int num_vertices = 2;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( num_vertices*my_rank + i );
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
    for ( int i = 0; i < num_vertices; ++i )
    {
	line_connectivity.push_back( vertex_handles[i] );
    }

    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> line_handle_array( line_handles.size() );
    std::copy( line_handles.begin(), line_handles.end(), 
	       line_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( line_connectivity.size() );
    std::copy( line_connectivity.begin(), line_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_LINE_SEGMENT, num_vertices,
				line_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Tri mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildTriContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 2;
    int num_vertices = 3;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( num_vertices*my_rank + i );
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
    for ( int i = 0; i < num_vertices; ++i )
    {
	tri_connectivity.push_back( vertex_handles[i] );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> tri_handle_array( tri_handles.size() );
    std::copy( tri_handles.begin(), tri_handles.end(), 
	       tri_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( tri_connectivity.size() );
    std::copy( tri_connectivity.begin(), tri_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_TRIANGLE, num_vertices,
				tri_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Quad mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildQuadContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 2;
    int num_vertices = 4;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( num_vertices*my_rank + i );
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
    for ( int i = 0; i < num_vertices; ++i )
    {
	quad_connectivity.push_back( vertex_handles[i] );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> quad_handle_array( quad_handles.size() );
    std::copy( quad_handles.begin(), quad_handles.end(), 
	       quad_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( quad_connectivity.size() );
    std::copy( quad_connectivity.begin(), quad_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_QUADRILATERAL, num_vertices,
				quad_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Shifted quad mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildShiftedQuadContainer( int my_rank,
							       int my_size )
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 2;
    int num_vertices = 4;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( num_vertices*my_rank + i + 3*my_size );
    }

    // x
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank );

    // y
    coords.push_back( my_rank-1 ); 
    coords.push_back( my_rank-1 ); 
    coords.push_back( my_rank ); 
    coords.push_back( my_rank ); 

    // Make the quadahedron.
    Teuchos::Array<int> quad_handles;
    Teuchos::Array<int> quad_connectivity;
    
    // handles
    quad_handles.push_back( 12 + my_size );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	quad_connectivity.push_back( vertex_handles[i] );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> quad_handle_array( quad_handles.size() );
    std::copy( quad_handles.begin(), quad_handles.end(), 
	       quad_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( quad_connectivity.size() );
    std::copy( quad_connectivity.begin(), quad_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_QUADRILATERAL, num_vertices,
				quad_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Tet mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildTetContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 3;
    int num_vertices = 4;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( num_vertices*my_rank + i );
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
    for ( int i = 0; i < num_vertices; ++i )
    {
	tet_connectivity.push_back( vertex_handles[i] );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> tet_handle_array( tet_handles.size() );
    std::copy( tet_handles.begin(), tet_handles.end(), 
	       tet_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( tet_connectivity.size() );
    std::copy( tet_connectivity.begin(), tet_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_TETRAHEDRON, num_vertices,
				tet_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Hex mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildHexContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 3;
    int num_vertices = 8;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( num_vertices*my_rank + i );
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
    for ( int i = 0; i < num_vertices; ++i )
    {
	hex_connectivity.push_back( vertex_handles[i] );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> hex_handle_array( hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(), 
	       hex_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_HEXAHEDRON, num_vertices,
				hex_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Shifted hex mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildShiftedHexContainer( int my_rank,
							      int my_size )
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 3;
    int num_vertices = 8;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( num_vertices*my_rank + i + my_size*4 );
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
    coords.push_back( my_rank-1 );
    coords.push_back( my_rank-1 );
    coords.push_back( my_rank-1 );
    coords.push_back( my_rank-1 );
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank );

    // Make the hexahedron.
    Teuchos::Array<int> hex_handles;
    Teuchos::Array<int> hex_connectivity;
    
    // handles
    hex_handles.push_back( 12+my_rank+my_size );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	hex_connectivity.push_back( vertex_handles[i] );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> hex_handle_array( hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(), 
	       hex_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_HEXAHEDRON, num_vertices,
				hex_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Pyramid mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildPyramidContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 3;
    int num_vertices = 5;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( num_vertices*my_rank + i );
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

    // Make the pyramid.
    Teuchos::Array<int> pyramid_handles;
    Teuchos::Array<int> pyramid_connectivity;
    
    // handles
    pyramid_handles.push_back( 12+my_rank );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	pyramid_connectivity.push_back( vertex_handles[i] );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> pyramid_handle_array( pyramid_handles.size() );
    std::copy( pyramid_handles.begin(), pyramid_handles.end(), 
	       pyramid_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( pyramid_connectivity.size() );
    std::copy( pyramid_connectivity.begin(), pyramid_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_PYRAMID, num_vertices,
				pyramid_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Shifted pyramid mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildShiftedPyramidContainer( int my_rank,
								  int my_size )
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 3;
    int num_vertices = 5;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( num_vertices*my_rank + i + 4*my_size + 8*my_size );
    }

    // x
    coords.push_back( my_rank ); 
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank+1 );
    coords.push_back( my_rank );

    // y
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank ); 
    coords.push_back( my_rank ); 

    // z
    coords.push_back( my_rank-1 );
    coords.push_back( my_rank-1 );
    coords.push_back( my_rank-1 );
    coords.push_back( my_rank-1 );
    coords.push_back( my_rank-2 );

    // Make the pyramid.
    Teuchos::Array<int> pyramid_handles;
    Teuchos::Array<int> pyramid_connectivity;
    
    // handles
    pyramid_handles.push_back( 12+2*my_size+my_rank );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	pyramid_connectivity.push_back( vertex_handles[i] );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> pyramid_handle_array( pyramid_handles.size() );
    std::copy( pyramid_handles.begin(), pyramid_handles.end(), 
	       pyramid_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( pyramid_connectivity.size() );
    std::copy( pyramid_connectivity.begin(), pyramid_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_PYRAMID, num_vertices,
				pyramid_handle_array, connectivity_array,
				permutation_list ) );
}

//---------------------------------------------------------------------------//
// Wedge mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildWedgeContainer( int my_rank )
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 3;
    int num_vertices = 6;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( i );
    }

    // x
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank );
    coords.push_back( my_rank+1 );
    coords.push_back( my_rank+1 );

    // y
    coords.push_back( my_rank ); 
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 
    coords.push_back( my_rank ); 
    coords.push_back( my_rank ); 
    coords.push_back( my_rank+1 ); 

    // z
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank );
    coords.push_back( my_rank+1 );
    coords.push_back( my_rank+1 );
    coords.push_back( my_rank+1 ); 

    // Make the wedge.
    Teuchos::Array<int> wedge_handles;
    Teuchos::Array<int> wedge_connectivity;
    
    // handles
    wedge_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	wedge_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<int> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<int> wedge_handle_array( wedge_handles.size() );
    std::copy( wedge_handles.begin(), wedge_handles.end(), 
	       wedge_handle_array.begin() );

    Teuchos::ArrayRCP<int> connectivity_array( wedge_connectivity.size() );
    std::copy( wedge_connectivity.begin(), wedge_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer<int>( vertex_dim, vertex_handle_array, coords_array,
				DTK_WEDGE, num_vertices,
				wedge_handle_array, connectivity_array,
				permutation_list ) );
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildLineContainer( my_rank );

    // Create a mesh manager.
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 1 ) );

    // Create a rendezvous.
    Rendezvous<MeshType> rendezvous( getDefaultComm<int>(), mesh_manager->dim(), box );
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
    Teuchos::Array<int> destinations = rendezvous.procsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points and check
    // that they are found in the correct element.
    int num_found = 0;
    Teuchos::Array<int> elements, elem_src_procs;
    rendezvous.elementsContainingPoints( points, elements, elem_src_procs );

    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] )
	{
	    TEST_ASSERT( elements[i] == std::numeric_limits<int>::max() );
	}
	else if ( elements[i] != std::numeric_limits<int>::max() )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) );
	    ++num_found;
	}
    }
    TEST_ASSERT( num_found > 0 );
}

//---------------------------------------------------------------------------//
// Tri mesh.
TEUCHOS_UNIT_TEST( MeshContainer, tri_rendezvous_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Create a bounding box that covers the entire mesh.
    double min = -Teuchos::ScalarTraits<double>::rmax();
    double max = Teuchos::ScalarTraits<double>::rmax();
    BoundingBox box( -100, -100, min, 100, 100, max );

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildTriContainer( my_rank );

    // Create a mesh manager. 
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 2 ) );

    // Create a rendezvous.
    Rendezvous<MeshType> rendezvous( getDefaultComm<int>(), mesh_manager->dim(), box );
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
    Teuchos::Array<int> destinations = rendezvous.procsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    int num_found = 0;
    Teuchos::Array<int> elements, elem_src_procs;
    rendezvous.elementsContainingPoints( points, elements, elem_src_procs );

    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] ||
	     points[num_points + i] < 0.0 || my_size < points[num_points + i] )
	{
	    TEST_ASSERT( elements[i] == std::numeric_limits<int>::max() );
	}
	else if ( elements[i] != std::numeric_limits<int>::max() )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) );
	    ++num_found;
	}
    }
    TEST_ASSERT( num_found > 0 );
}

//---------------------------------------------------------------------------//
// Quad mesh.
TEUCHOS_UNIT_TEST( MeshContainer, quad_rendezvous_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Create a bounding box that covers the entire mesh.
    double min = -Teuchos::ScalarTraits<double>::rmax();
    double max = Teuchos::ScalarTraits<double>::rmax();
    BoundingBox box( -100, -100, min, 100, 100, max );

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildQuadContainer( my_rank );

    // Create a mesh manager. 
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 2 ) );

    // Create a rendezvous.
    Rendezvous<MeshType> rendezvous( getDefaultComm<int>(), mesh_manager->dim(), box );
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
    Teuchos::Array<int> destinations = rendezvous.procsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    int num_found = 0;
    Teuchos::Array<int> elements, elem_src_procs;
    rendezvous.elementsContainingPoints( points, elements, elem_src_procs );

    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] ||
	     points[num_points + i] < 0.0 || my_size < points[num_points + i] )
	{
	    TEST_ASSERT( elements[i] == std::numeric_limits<int>::max() );
	}
	else if ( elements[i] != std::numeric_limits<int>::max() )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) );
	    ++num_found;
	}
    }
    TEST_ASSERT( num_found > 0 );
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildTetContainer( my_rank );

    // Create a mesh manager. 
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 3 ) );

    // Create a rendezvous.
    Rendezvous<MeshType> rendezvous( getDefaultComm<int>(), mesh_manager->dim(), box );
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
    Teuchos::Array<int> destinations = rendezvous.procsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    Teuchos::Array<int> elements, elem_src_procs;
    rendezvous.elementsContainingPoints( points, elements, elem_src_procs );

    int num_found = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] ||
	     points[num_points + i] < 0.0 || my_size < points[num_points + i] ||
	     points[2*num_points + i] < 0.0 || my_size < points[2*num_points + i] )
	{
	    TEST_ASSERT( elements[i] == std::numeric_limits<int>::max() );
	}
	else if ( elements[i] != std::numeric_limits<int>::max() )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) );
	    ++num_found;
	}
    }
    TEST_ASSERT( num_found > 0 );
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildHexContainer( my_rank );

    // Create a mesh manager. 
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 3 ) );

    // Create a rendezvous.
    Rendezvous<MeshType> rendezvous( getDefaultComm<int>(), mesh_manager->dim(), box );
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
    Teuchos::Array<int> destinations = rendezvous.procsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    Teuchos::Array<int> elements, elem_src_procs;
    rendezvous.elementsContainingPoints( points, elements, elem_src_procs );

    int num_found = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] ||
	     points[num_points + i] < 0.0 || my_size < points[num_points + i] ||
	     points[2*num_points + i] < 0.0 || my_size < points[2*num_points + i] )
	{
	    TEST_ASSERT( elements[i] == std::numeric_limits<int>::max() );
	}
	else if ( elements[i] != std::numeric_limits<int>::max() )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) );
	    ++num_found;
	}
    }
    TEST_ASSERT( num_found > 0 );
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildPyramidContainer( my_rank );

    // Create a mesh manager. 
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 3 ) );

    // Create a rendezvous.
    Rendezvous<MeshType> rendezvous( getDefaultComm<int>(), mesh_manager->dim(), box );
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
    Teuchos::Array<int> destinations = rendezvous.procsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    Teuchos::Array<int> elements, elem_src_procs;
    rendezvous.elementsContainingPoints( points, elements, elem_src_procs );

    int num_found = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] ||
	     points[num_points + i] < 0.0 || my_size < points[num_points + i] ||
	     points[2*num_points + i] < 0.0 || my_size < points[2*num_points + i] )
	{
	    TEST_ASSERT( elements[i] == std::numeric_limits<int>::max() );
	}
	else if ( elements[i] != std::numeric_limits<int>::max() )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) );
	    ++num_found;
	}
    }
    TEST_ASSERT( num_found > 0 );
}

//---------------------------------------------------------------------------//
// Wedge mesh.
TEUCHOS_UNIT_TEST( MeshContainer, wedge_rendezvous_test )
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildWedgeContainer( my_rank );

    // Create a mesh manager. 
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 3 ) );

    // Create a rendezvous.
    Rendezvous<MeshType> rendezvous( getDefaultComm<int>(), mesh_manager->dim(), box );
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
    Teuchos::Array<int> destinations = rendezvous.procsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    int num_found = 0;
    Teuchos::Array<int> elements, elem_src_procs;
    rendezvous.elementsContainingPoints( points, elements, elem_src_procs );

    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] ||
	     points[num_points + i] < 0.0 || my_size < points[num_points + i] ||
	     points[2*num_points+i] < 0.0 || points[2*num_points+i] > 1.0 )
	{
	    TEST_ASSERT( elements[i] == std::numeric_limits<int>::max() );
	}
	else if ( elements[i] != std::numeric_limits<int>::max() )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) );
	    ++num_found;
	}
    }
    TEST_ASSERT( num_found > 0 );
}

//---------------------------------------------------------------------------//
// 2d hybrid test.
TEUCHOS_UNIT_TEST( MeshContainer, 2d_hybrid_rendezvous_test )
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 2 );
    mesh_blocks[0] = buildTriContainer( my_rank );
    mesh_blocks[1] = buildShiftedQuadContainer( my_rank, my_size );

    // Create a mesh manager. 
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 2 ) );

    // Create a rendezvous.
    Rendezvous<MeshType> rendezvous( getDefaultComm<int>(), mesh_manager->dim(), box );
    rendezvous.build( mesh_manager );

    // Give every process a unique block of random numbers;
    std::srand( my_rank*num_points*mesh_manager->dim() );

    // Create some random points.
    int num_rand = num_points*mesh_manager->dim();
    Teuchos::ArrayRCP<double> points( num_rand );
    for ( int i = 0; i < num_rand; ++i )
    {
	points[i] = (my_size+2) * (double) std::rand() / RAND_MAX - 0.5;
    }

    // Get the destination procs for the random points.
    Teuchos::Array<int> destinations = rendezvous.procsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    Teuchos::Array<int> elements, elem_src_procs;
    rendezvous.elementsContainingPoints( points, elements, elem_src_procs );

    int num_found = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] ||
	     points[num_points + i] < -1.0 || my_size < points[num_points + i] )
	{
	    TEST_ASSERT( elements[i] == std::numeric_limits<int>::max() );
	}
	else if ( elements[i] != std::numeric_limits<int>::max() )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) ||
			 elements[i] == 12 + std::floor( points[i] ) + my_size );
	    ++num_found;
	}
    }
}

//---------------------------------------------------------------------------//
// 3d hybrid test.
TEUCHOS_UNIT_TEST( MeshContainer, 3d_hybrid_rendezvous_test )
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 3 );
    mesh_blocks[0] = buildTetContainer( my_rank );
    mesh_blocks[1] = buildShiftedPyramidContainer( my_rank, my_size );
    mesh_blocks[2] = buildShiftedHexContainer( my_rank, my_size );

    // Create a mesh manager. 
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp(
	new MeshManager<MeshType>( mesh_blocks, getDefaultComm<int>(), 3 ) );

    // Create a rendezvous.
    Rendezvous<MeshType> rendezvous( getDefaultComm<int>(), mesh_manager->dim(), box );
    rendezvous.build( mesh_manager );

    // Give every process a unique block of random numbers;
    std::srand( my_rank*num_points*mesh_manager->dim() );

    // Create some random points.
    int num_rand = num_points*mesh_manager->dim();
    Teuchos::ArrayRCP<double> points( num_rand );
    for ( int i = 0; i < num_rand; ++i )
    {
	points[i] = (my_size+3) * (double) std::rand() / RAND_MAX - 2.5;
    }

    // Get the destination procs for the random points.
    Teuchos::Array<int> destinations = rendezvous.procsContainingPoints( points );
    for ( int i = 0; i < num_points; ++i )
    {
	TEST_ASSERT( 0 <= destinations[i] && destinations[i] < my_size );
    }

    // Search the rendezvous decomposition for some random points.
    Teuchos::Array<int> elements, elem_src_procs;
    rendezvous.elementsContainingPoints( points, elements, elem_src_procs );

    int num_found = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	if ( points[i] < 0.0 || my_size < points[i] ||
	     points[num_points + i] < 0.0 || my_size < points[num_points + i] ||
	     points[2*num_points + i] < -2.0 || my_size < points[2*num_points + i] )
	{
	    TEST_ASSERT( elements[i] == std::numeric_limits<int>::max() );
	}
	else if ( elements[i] != std::numeric_limits<int>::max() )
	{
	    TEST_ASSERT( elements[i] == 12 + std::floor( points[i] ) ||
			 elements[i] == 12 + std::floor( points[i] ) + my_size ||
			 elements[i] == 12 + std::floor( points[i] ) + 2*my_size );
	    ++num_found;
	}
    }
    TEST_ASSERT( num_found > 0 );
}

//---------------------------------------------------------------------------//
// end tstRendezvous.cpp
//---------------------------------------------------------------------------//
