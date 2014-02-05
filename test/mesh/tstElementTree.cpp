//---------------------------------------------------------------------------//
/*! 
 * \file tstElementTree.cpp
 * \author Stuart R. Slattery
 * \brief ElementTree unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_ElementTree.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshManager.hpp>
#include <DTK_MeshTools.hpp>
#include <DTK_MeshContainer.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>

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
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildLineContainer()
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
	vertex_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 

    // Make the line.
    Teuchos::Array<int> line_handles;
    Teuchos::Array<int> line_connectivity;
    
    // handles
    line_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	line_connectivity.push_back( i );
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
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildTriContainer()
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
	vertex_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 

    // Make the triahedron.
    Teuchos::Array<int> tri_handles;
    Teuchos::Array<int> tri_connectivity;
    
    // handles
    tri_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	tri_connectivity.push_back( i );
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
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildQuadContainer()
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
	vertex_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 );

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 

    // Make the quadahedron.
    Teuchos::Array<int> quad_handles;
    Teuchos::Array<int> quad_connectivity;
    
    // handles
    quad_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	quad_connectivity.push_back( i );
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
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildShiftedQuadContainer()
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
	vertex_handles.push_back( i );
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
    for ( int i = 0; i < num_vertices; ++i )
    {
	quad_connectivity.push_back( i );
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
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildTetContainer()
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
	vertex_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 );

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 ); 

    // z
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );

    // Make the tetahedron.
    Teuchos::Array<int> tet_handles;
    Teuchos::Array<int> tet_connectivity;
    
    // handles
    tet_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	tet_connectivity.push_back( i );
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
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildHexContainer()
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
	vertex_handles.push_back( i );
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
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );

    // Make the hexahedron.
    Teuchos::Array<int> hex_handles;
    Teuchos::Array<int> hex_connectivity;
    
    // handles
    hex_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	hex_connectivity.push_back( i );
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
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildShiftedHexContainer()
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
	vertex_handles.push_back( i );
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
    for ( int i = 0; i < num_vertices; ++i )
    {
	hex_connectivity.push_back( i );
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
// Parallel hex mesh.
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildParallelHexContainer()
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();

    // Make some vertices.
    Teuchos::Array<int> vertex_handles;
    Teuchos::Array<double> coords;

    int vertex_dim = 3;
    int num_vertices = 8;

    // handles
    for ( int i = 0; i < num_vertices; ++i )
    {
	vertex_handles.push_back( i );
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
    hex_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	hex_connectivity.push_back( i );
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
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildPyramidContainer()
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
	vertex_handles.push_back( i );
    }

    // x
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 ); 

    // z
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );

    // Make the pyramid.
    Teuchos::Array<int> pyramid_handles;
    Teuchos::Array<int> pyramid_connectivity;
    
    // handles
    pyramid_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	pyramid_connectivity.push_back( i );
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
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildShiftedPyramidContainer()
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
	vertex_handles.push_back( i );
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
    for ( int i = 0; i < num_vertices; ++i )
    {
	pyramid_connectivity.push_back( i );
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
Teuchos::RCP<DataTransferKit::MeshContainer<int> > buildWedgeContainer()
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
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 

    // z
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 );
    coords.push_back( 1.0 ); 

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
// Tests
//---------------------------------------------------------------------------//
// Line mesh.
TEUCHOS_UNIT_TEST( MeshContainer, line_element_tree_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildLineContainer();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>(mesh_blocks, getDefaultComm<int>(), 1) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree<MeshType> element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(1);
    int ordinal = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	ordinal = 0;
	point[0] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	if ( 0.0 <= point[0] && point[0] <= 1.0 )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 12 );
	}
	else
	{
	    TEST_ASSERT( !element_tree.findPoint( point, ordinal ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Tri mesh.
TEUCHOS_UNIT_TEST( MeshContainer, tri_element_tree_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildTriContainer();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>(mesh_blocks, getDefaultComm<int>(), 2) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree<MeshType> element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(2);
    int ordinal = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	ordinal = 0;
	point[0] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[1] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;

	if ( 0.0 <= point[0] && point[0] <= 1.0 &&
	     0.0 <= point[1] && point[1] <= point[0] )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 12 );
	}
	else
	{
	    TEST_ASSERT( !element_tree.findPoint( point, ordinal ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Quad mesh.
TEUCHOS_UNIT_TEST( MeshContainer, quad_element_tree_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildQuadContainer();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>(mesh_blocks, getDefaultComm<int>(), 2) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree<MeshType> element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(2);
    int ordinal = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	ordinal = 0;
	point[0] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[1] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;

	if ( 0.0 <= point[0] && point[0] <= 1.0 &&
	     0.0 <= point[1] && point[1] <= 1.0 )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 12 );
	}
	else
	{
	    TEST_ASSERT( !element_tree.findPoint( point, ordinal ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Tet mesh.
TEUCHOS_UNIT_TEST( MeshContainer, tet_element_tree_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildTetContainer();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>(mesh_blocks, getDefaultComm<int>(), 3) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree<MeshType> element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(3);
    int ordinal = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	ordinal = 0;
	point[0] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[1] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[2] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
    	
	if ( std::max( std::max(-point[0],-point[1]),
		       std::max(-point[2], point[0]+point[1]+point[2]-1) )
	     < 1.0e-8 )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 12 );
	}
	else
	{
	    TEST_ASSERT( !element_tree.findPoint( point, ordinal ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Hex mesh.
TEUCHOS_UNIT_TEST( MeshContainer, hex_element_tree_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildHexContainer();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>(mesh_blocks, getDefaultComm<int>(), 3) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree<MeshType> element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(3);
    int ordinal = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	ordinal = 0;
	point[0] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[1] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[2] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;

	if ( 0.0 <= point[0] && point[0] <= 1.0 &&
	     0.0 <= point[1] && point[1] <= 1.0 &&
	     0.0 <= point[2] && point[2] <= 1.0 )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 12 );
	}
	else
	{
	    TEST_ASSERT( !element_tree.findPoint( point, ordinal ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Pyramid mesh.
TEUCHOS_UNIT_TEST( MeshContainer, pyramid_element_tree_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildPyramidContainer();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>(mesh_blocks, getDefaultComm<int>(), 3) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree<MeshType> element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(3);
    int ordinal = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	ordinal = 0;
	point[0] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[1] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[2] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
    	
	if( 0.0 <= point[0] && point[0] <= 1.0-point[2] &&
	    0.0 <= point[1] && point[1] <= 1.0-point[2] && 
	    0.0 <= point[2] && point[2] <= 1.0 )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 12 );
	}
	else
	{
	    TEST_ASSERT( !element_tree.findPoint( point, ordinal ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Wedge mesh.
TEUCHOS_UNIT_TEST( MeshContainer, wedge_element_tree_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildWedgeContainer();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>(mesh_blocks, getDefaultComm<int>(), 3) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree<MeshType> element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(3);
    int ordinal = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	ordinal = 0;
	point[0] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[1] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[2] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;

	if ( 0.0 <= point[0] && point[0] <= 1.0 &&
	     0.0 <= point[1] && point[1] <= point[0] &&
	     0.0 <= point[2] && point[2] <= 1.0 )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 12 );
	}
	else
	{
	    TEST_ASSERT( !element_tree.findPoint( point, ordinal ) );
	}
    }
}

//---------------------------------------------------------------------------//
// Parallel hex mesh.
TEUCHOS_UNIT_TEST( MeshContainer, parallel_hex_element_tree_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = buildParallelHexContainer();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>(mesh_blocks, getDefaultComm<int>(), 3) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree<MeshType> element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(3);

    int ordinal = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	ordinal = 0;
	point[0] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[1] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[2] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;

	if ( 0.0 <= point[0] && point[0] <= 1.0 &&
	     0.0 <= point[1] && point[1] <= 1.0 &&
	     my_rank <= point[2] && point[2] <= my_rank+1 )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 12 );
	}
	else
	{
	    TEST_ASSERT( !element_tree.findPoint( point, ordinal ) );
	}
    }
}

//---------------------------------------------------------------------------//
// 2d hybrid test.
TEUCHOS_UNIT_TEST( MeshContainer, 2d_hybrid_element_tree_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 2 );
    mesh_blocks[0] = buildTriContainer();
    mesh_blocks[1] = buildShiftedQuadContainer();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>(mesh_blocks, getDefaultComm<int>(), 2) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree<MeshType> element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(2);
    int ordinal = 0;
    double tol = 1.0e-8;
    for ( int i = 0; i < num_points; ++i )
    {
	ordinal = 0;
	point[0] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[1] = 3.0 * (double) std::rand() / RAND_MAX - 1.5;

	// We can end up either in the quad or tri on their boundary.
	if ( 0.0 <= point[0] && point[0] <= 1.0 &&
	     -tol <= point[1] && point[1] <= tol )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 12 || ordinal == 9 );
	}
	// In the tri.
	else if ( 0.0 <= point[0] && point[0] <= 1.0 &&
		  0.0 < point[1] && point[1] <= point[0] )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 12 );
	}
	// In the quad.
	else if ( 0.0 <= point[0] && point[0] <= 1.0 &&
		  -1.0 <= point[1] && point[1] < 0.0 )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 9 );
	}
	// Neither
	else
	{
	    TEST_ASSERT( !element_tree.findPoint( point, ordinal ) );
	}
    }
}

//---------------------------------------------------------------------------//
// 3d hybrid test.
TEUCHOS_UNIT_TEST( MeshContainer, 3d_hybrid_element_tree_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 3 );
    mesh_blocks[0] = buildTetContainer();
    mesh_blocks[1] = buildShiftedPyramidContainer();
    mesh_blocks[2] = buildShiftedHexContainer();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>(mesh_blocks, getDefaultComm<int>(), 3) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree<MeshType> element_tree( mesh_manager );

    // Search the tree for some random points.
    double tol = 1.0e-8;
    int num_points = 1000;
    Teuchos::Array<double> point(3);
    int ordinal = 0;
    for ( int i = 0; i < num_points; ++i )
    {
	ordinal = 0;
	point[0] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[1] = 2.0 * (double) std::rand() / RAND_MAX - 0.5;
	point[2] = 4.0 * (double) std::rand() / RAND_MAX - 2.5;
    	
	// Hex/Pyramid boundary.
	if ( 0.0 <= point[0] && point[0] <= 1.0 &&
	     0.0 <= point[1] && point[1] <= 1.0 &&
	     -1.0-tol <= point[2] && point[2] <= -1.0+tol )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 6 || ordinal == 89 );
	}

	// Hex/Tet boundary.
	else if ( 0.0 <= point[0] && point[0] <= 1.0 &&
		  0.0 <= point[1] && point[1] <= 1.0 &&
		  -tol <= point[2] && point[2] <= tol )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 12 || ordinal == 6 );
	}

	// Tet
	else if ( std::max( std::max(-point[0],-point[1]),
			    std::max(-point[2], point[0]+point[1]+point[2]-1) )
		  < tol )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 12 );
	}

	// Hex
	else if ( 0.0 <= point[0] && point[0] <= 1.0 &&
		  0.0 <= point[1] && point[1] <= 1.0 && 
		  -1.0 <= point[2] && point[2] <= 0.0 )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 6 );
	}
	// Pyramid
	else if( 0.0 <= point[0] && point[0] <= 2.0+point[2] &&
		 0.0 <= point[1] && point[1] <= 2.0+point[2] && 
		 -2.0 <= point[2] && point[2] <= -1.0 )
	{
	    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
	    TEST_ASSERT( ordinal == 89 );
	}

	// None
	else
	{
	    TEST_ASSERT( !element_tree.findPoint( point, ordinal ) );
	}
    }
}

//---------------------------------------------------------------------------//
// end tstElementTree.cpp
//---------------------------------------------------------------------------//
