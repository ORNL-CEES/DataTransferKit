//---------------------------------------------------------------------------//
/*!
 * \file MeshBuilder.hpp
 * \author Stuart R. Slattery
 * \brief Meshes for unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_MeshContainer.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshBlock.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_TypeTraits.hpp>

namespace DataTransferKit
{
namespace UnitTest
{
namespace MeshBuilder
{
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
Teuchos::RCP<DataTransferKit::MeshBlock> buildLineBlock()
{
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
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array(
	vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> line_handle_array( 
	line_handles.size() );
    std::copy( line_handles.begin(), line_handles.end(), 
	       line_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array( 
	line_connectivity.size() );
    std::copy( line_connectivity.begin(), line_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( new DataTransferKit::MeshContainer( 
			     vertex_dim, vertex_handle_array, coords_array,
			     DataTransferKit::DTK_LINE_SEGMENT, num_vertices,
			     line_handle_array, connectivity_array,
			     permutation_list) );
}

//---------------------------------------------------------------------------//
// Tri mesh.
Teuchos::RCP<DataTransferKit::MeshBlock> buildTriBlock()
{
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
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array( 
	vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> tri_handle_array( 
	tri_handles.size() );
    std::copy( tri_handles.begin(), tri_handles.end(), 
	       tri_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array( 
	tri_connectivity.size() );
    std::copy( tri_connectivity.begin(), tri_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( new DataTransferKit::MeshContainer( 
			     vertex_dim, vertex_handle_array, coords_array,
			     DataTransferKit::DTK_TRIANGLE, num_vertices,
			     tri_handle_array, connectivity_array,
			     permutation_list)  );
}

//---------------------------------------------------------------------------//
// Quad mesh.
Teuchos::RCP<DataTransferKit::MeshBlock> buildQuadBlock()
{
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
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array( 
	vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> quad_handle_array(
	quad_handles.size() );
    std::copy( quad_handles.begin(), quad_handles.end(), 
	       quad_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array( 
	quad_connectivity.size() );
    std::copy( quad_connectivity.begin(), quad_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( new DataTransferKit::MeshContainer( 
			     vertex_dim, vertex_handle_array, coords_array,
			     DataTransferKit::DTK_QUADRILATERAL, num_vertices,
			     quad_handle_array, connectivity_array,
			     permutation_list) );
}

//---------------------------------------------------------------------------//
// Tet mesh.
Teuchos::RCP<DataTransferKit::MeshBlock> buildTetBlock()
{
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

    // Make the tetrahedron.
    Teuchos::Array<int> tet_handles;
    Teuchos::Array<int> tet_connectivity;
    
    // handles
    tet_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	tet_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array( 
	vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> tet_handle_array( 
	tet_handles.size() );
    std::copy( tet_handles.begin(), tet_handles.end(), 
	       tet_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array( 
	tet_connectivity.size() );
    std::copy( tet_connectivity.begin(), tet_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( new DataTransferKit::MeshContainer( 
			     vertex_dim, vertex_handle_array, coords_array,
			     DataTransferKit::DTK_TETRAHEDRON, num_vertices,
			     tet_handle_array, connectivity_array,
			     permutation_list) );
}

//---------------------------------------------------------------------------//
// Hex mesh.
Teuchos::RCP<DataTransferKit::MeshBlock> buildHexBlock()
{
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
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array(
	vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> hex_handle_array( 
	hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(), 
	       hex_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array( 
	hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( new DataTransferKit::MeshContainer( 
			     vertex_dim, vertex_handle_array, coords_array,
			     DataTransferKit::DTK_HEXAHEDRON, num_vertices,
			     hex_handle_array, connectivity_array,
			     permutation_list) );
}

//---------------------------------------------------------------------------//
// Pyramid mesh.
Teuchos::RCP<DataTransferKit::MeshBlock> buildPyramidBlock()
{
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

    // Make the pyramid
    Teuchos::Array<int> pyramid_handles;
    Teuchos::Array<int> pyramid_connectivity;
    
    // handles
    pyramid_handles.push_back( 12 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	pyramid_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array(
	vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> pyramid_handle_array( 
	pyramid_handles.size() );
    std::copy( pyramid_handles.begin(), pyramid_handles.end(), 
	       pyramid_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array(
	pyramid_connectivity.size() );
    std::copy( pyramid_connectivity.begin(), pyramid_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( new DataTransferKit::MeshContainer( 
			     vertex_dim, vertex_handle_array, coords_array,
			     DataTransferKit::DTK_PYRAMID, num_vertices,
			     pyramid_handle_array, connectivity_array,
			     permutation_list) );
}

//---------------------------------------------------------------------------//
// Wedge mesh.
Teuchos::RCP<DataTransferKit::MeshBlock> buildWedgeBlock()
{
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
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array( 
	vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> wedge_handle_array( 
	wedge_handles.size() );
    std::copy( wedge_handles.begin(), wedge_handles.end(), 
	       wedge_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array(
	wedge_connectivity.size() );
    std::copy( wedge_connectivity.begin(), wedge_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( new DataTransferKit::MeshContainer( 
			     vertex_dim, vertex_handle_array, coords_array,
			     DataTransferKit::DTK_WEDGE, num_vertices,
			     wedge_handle_array, connectivity_array,
			     permutation_list) );
}


//---------------------------------------------------------------------------//
// Shifted quad mesh.
Teuchos::RCP<DataTransferKit::MeshBlock> buildShiftedQuadBlock()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<MeshId> vertex_handles;
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
    Teuchos::Array<MeshId> quad_handles;
    Teuchos::Array<MeshId> quad_connectivity;
    
    // handles
    quad_handles.push_back( 9 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	quad_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<MeshId> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<MeshId> quad_handle_array( quad_handles.size() );
    std::copy( quad_handles.begin(), quad_handles.end(), 
	       quad_handle_array.begin() );

    Teuchos::ArrayRCP<MeshId> connectivity_array( quad_connectivity.size() );
    std::copy( quad_connectivity.begin(), quad_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer( vertex_dim, vertex_handle_array, coords_array,
			   DTK_QUADRILATERAL, num_vertices,
			   quad_handle_array, connectivity_array,
			   permutation_list ) );
}
//---------------------------------------------------------------------------//
// Parallel hex mesh.
Teuchos::RCP<DataTransferKit::MeshBlock> buildParallelHexBlock()
{
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
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array( 
	vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> hex_handle_array( 
	hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(), 
	       hex_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array( 
	hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( new DataTransferKit::MeshContainer(
			     vertex_dim, vertex_handle_array, coords_array,
			     DataTransferKit::DTK_HEXAHEDRON, num_vertices,
			     hex_handle_array, connectivity_array,
			     permutation_list) );
}

//---------------------------------------------------------------------------//
// Shifted hex mesh.
Teuchos::RCP<DataTransferKit::MeshBlock> buildShiftedHexBlock()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<MeshId> vertex_handles;
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
    Teuchos::Array<MeshId> hex_handles;
    Teuchos::Array<MeshId> hex_connectivity;
    
    // handles
    hex_handles.push_back( 6 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	hex_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<MeshId> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<MeshId> hex_handle_array( hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(), 
	       hex_handle_array.begin() );

    Teuchos::ArrayRCP<MeshId> connectivity_array( hex_connectivity.size() );
    std::copy( hex_connectivity.begin(), hex_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp(
	new MeshContainer( vertex_dim, vertex_handle_array, coords_array,
			   DTK_HEXAHEDRON, num_vertices,
			   hex_handle_array, connectivity_array,
			   permutation_list ) );
}

//---------------------------------------------------------------------------//
// Shifted pyramid mesh.
Teuchos::RCP<DataTransferKit::MeshBlock> buildShiftedPyramidBlock()
{
    using namespace DataTransferKit;

    // Make some vertices.
    Teuchos::Array<MeshId> vertex_handles;
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
    Teuchos::Array<MeshId> pyramid_handles;
    Teuchos::Array<MeshId> pyramid_connectivity;
    
    // handles
    pyramid_handles.push_back( 89 );

    // connectivity
    for ( int i = 0; i < num_vertices; ++i )
    {
	pyramid_connectivity.push_back( i );
    }
    
    Teuchos::ArrayRCP<MeshId> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<MeshId> pyramid_handle_array( pyramid_handles.size() );
    std::copy( pyramid_handles.begin(), pyramid_handles.end(), 
	       pyramid_handle_array.begin() );

    Teuchos::ArrayRCP<MeshId> connectivity_array( pyramid_connectivity.size() );
    std::copy( pyramid_connectivity.begin(), pyramid_connectivity.end(), 
	       connectivity_array.begin() );

    Teuchos::ArrayRCP<int> permutation_list( num_vertices );
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i;
    }
    
    return Teuchos::rcp( 
	new MeshContainer( vertex_dim, vertex_handle_array, coords_array,
			   DTK_PYRAMID, num_vertices,
			   pyramid_handle_array, connectivity_array,
			   permutation_list ) );
}

//---------------------------------------------------------------------------//

} // end namespace MeshBuilder
} // end namespace UnitTest
} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end MeshBuilder.hpp
//---------------------------------------------------------------------------//
