//---------------------------------------------------------------------------//
/*!
 * \file tstMeshContainer.cpp
 * \author Stuart R. Slattery
 * \brief MeshContainer unit tests.
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
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> line_handle_array( line_handles.size() );
    std::copy( line_handles.begin(), line_handles.end(), 
	       line_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array( line_connectivity.size() );
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
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> tri_handle_array( tri_handles.size() );
    std::copy( tri_handles.begin(), tri_handles.end(), 
	       tri_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array( tri_connectivity.size() );
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
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> quad_handle_array( quad_handles.size() );
    std::copy( quad_handles.begin(), quad_handles.end(), 
	       quad_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array( quad_connectivity.size() );
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
    coords.push_back( 1.0 ); 
    coords.push_back( 0.0 );

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 

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
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> tet_handle_array( tet_handles.size() );
    std::copy( tet_handles.begin(), tet_handles.end(), 
	       tet_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array( tet_connectivity.size() );
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
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> hex_handle_array( hex_handles.size() );
    std::copy( hex_handles.begin(), hex_handles.end(), 
	       hex_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array( hex_connectivity.size() );
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
    coords.push_back( 0.5 );

    // y
    coords.push_back( 0.0 ); 
    coords.push_back( 0.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 1.0 ); 
    coords.push_back( 0.5 ); 

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
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> pyramid_handle_array( pyramid_handles.size() );
    std::copy( pyramid_handles.begin(), pyramid_handles.end(), 
	       pyramid_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array( pyramid_connectivity.size() );
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
    coords.push_back( 0.5 ); 
    coords.push_back( 0.0 );
    coords.push_back( 1.0 );
    coords.push_back( 0.5 );

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
    
    Teuchos::ArrayRCP<DataTransferKit::MeshId> vertex_handle_array( vertex_handles.size() );
    std::copy( vertex_handles.begin(), vertex_handles.end(), 
	       vertex_handle_array.begin() );

    Teuchos::ArrayRCP<double> coords_array( coords.size() );
    std::copy( coords.begin(), coords.end(), coords_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> wedge_handle_array( wedge_handles.size() );
    std::copy( wedge_handles.begin(), wedge_handles.end(), 
	       wedge_handle_array.begin() );

    Teuchos::ArrayRCP<DataTransferKit::MeshId> connectivity_array( wedge_connectivity.size() );
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
// Tests
//---------------------------------------------------------------------------//
// Line mesh.
TEUCHOS_UNIT_TEST( MeshBlock, line_block_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = buildLineBlock();

    // Mesh parameters.
    int vertex_dim = 1;
    int num_vertices = 2;
    int element_topo = DataTransferKit::DTK_LINE_SEGMENT;

    // Basic block info.
    TEST_ASSERT( (int) mesh_block->dimension() == vertex_dim );
    TEST_ASSERT( (int) mesh_block->verticesPerElement() == num_vertices );
    TEST_ASSERT( (int) mesh_block->elementTopology() == element_topo );

    // Vertices.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->vertexIds()[i] == i );
    }

    // Coords.
    // x
    TEST_ASSERT( mesh_block->vertexCoordinates()[0] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[1] == 1.0 ); 

    // Elements.
    TEST_ASSERT( mesh_block->elementIds()[0] == 12 );

    // Connectivity.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->connectivity()[i] == i );
    }

    // Permutation.
    for ( int i = 0; i < num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->permutation()[i] == i );
    }
}

//---------------------------------------------------------------------------//
// Tri mesh.
TEUCHOS_UNIT_TEST( MeshBlock, tri_block_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = buildTriBlock();

    // Mesh parameters.
    int vertex_dim = 2;
    int num_vertices = 3;
    int element_topo = DataTransferKit::DTK_TRIANGLE;

    // Basic block info.
    TEST_ASSERT( (int) mesh_block->dimension() == vertex_dim );
    TEST_ASSERT( (int) mesh_block->verticesPerElement() == num_vertices );
    TEST_ASSERT( (int) mesh_block->elementTopology() == element_topo );

    // Vertices.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->vertexIds()[i] == i );
    }

    // Coords.
    // x
    TEST_ASSERT( mesh_block->vertexCoordinates()[0] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[1] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[2] == 1.0 ); 

    // y
    TEST_ASSERT( mesh_block->vertexCoordinates()[3] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[4] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[5] == 1.0 ); 

    // Elements.
    TEST_ASSERT( mesh_block->elementIds()[0] == 12 );

    // Connectivity.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->connectivity()[i] == i );
    }

    // Permutation.
    for ( int i = 0; i < num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->permutation()[i] == i );
    }
}

//---------------------------------------------------------------------------//
// Quad mesh.
TEUCHOS_UNIT_TEST( MeshBlock, quad_block_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = buildQuadBlock();

    // Mesh parameters.
    int vertex_dim = 2;
    int num_vertices = 4;
    int element_topo = DataTransferKit::DTK_QUADRILATERAL;

    // Basic block info.
    TEST_ASSERT( (int) mesh_block->dimension() == vertex_dim );
    TEST_ASSERT( (int) mesh_block->verticesPerElement() == num_vertices );
    TEST_ASSERT( (int) mesh_block->elementTopology() == element_topo );

    // Vertices.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->vertexIds()[i] == i );
    }

    // Coords.
    // x
    TEST_ASSERT( mesh_block->vertexCoordinates()[0] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[1] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[2] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[3] == 0.0 );

    // y
    TEST_ASSERT( mesh_block->vertexCoordinates()[4] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[5] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[6] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[7] == 1.0 ); 

    // Elements.
    TEST_ASSERT( mesh_block->elementIds()[0] == 12 );

    // Connectivity.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->connectivity()[i] == i );
    }

    // Permutation.
    for ( int i = 0; i < num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->permutation()[i] == i );
    }
}

//---------------------------------------------------------------------------//
// Tet mesh.
TEUCHOS_UNIT_TEST( MeshBlock, tet_block_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = buildTetBlock();

    // Mesh parameters.
    int vertex_dim = 3;
    int num_vertices = 4;
    int element_topo = DataTransferKit::DTK_TETRAHEDRON;

    // Basic block info.
    TEST_ASSERT( (int) mesh_block->dimension() == vertex_dim );
    TEST_ASSERT( (int) mesh_block->verticesPerElement() == num_vertices );
    TEST_ASSERT( (int) mesh_block->elementTopology() == element_topo );

    // Vertices.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->vertexIds()[i] == i );
    }

    // Coords.
    // x
    TEST_ASSERT( mesh_block->vertexCoordinates()[0] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[1] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[2] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[3] == 0.0 );

    // y
    TEST_ASSERT( mesh_block->vertexCoordinates()[4] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[5] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[6] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[7] == 1.0 ); 

    // z
    TEST_ASSERT( mesh_block->vertexCoordinates()[8]  == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[9]  == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[10] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[11] == 1.0 );

    // Elements.
    TEST_ASSERT( mesh_block->elementIds()[0] == 12 );

    // Connectivity.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->connectivity()[i] == i );
    }

    // Permutation.
    for ( int i = 0; i < num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->permutation()[i] == i );
    }
}

//---------------------------------------------------------------------------//
// Hex mesh.
TEUCHOS_UNIT_TEST( MeshBlock, hex_block_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = buildHexBlock();

    // Mesh parameters.
    int vertex_dim = 3;
    int num_vertices = 8;
    int element_topo = DataTransferKit::DTK_HEXAHEDRON;

    // Basic block info.
    TEST_ASSERT( (int) mesh_block->dimension() == vertex_dim );
    TEST_ASSERT( (int) mesh_block->verticesPerElement() == num_vertices );
    TEST_ASSERT( (int) mesh_block->elementTopology() == element_topo );

    // Vertices.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->vertexIds()[i] == i );
    }

    // Coords.
    // x
    TEST_ASSERT( mesh_block->vertexCoordinates()[0] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[1] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[2] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[3] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[4] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[5] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[6] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[7] == 0.0 ); 

    // y
    TEST_ASSERT( mesh_block->vertexCoordinates()[8]  == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[9]  == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[10] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[11] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[12] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[13] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[14] == 1.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[15] == 1.0 );

    // z
    TEST_ASSERT( mesh_block->vertexCoordinates()[16] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[17] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[18] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[19] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[20] == 1.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[21] == 1.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[22] == 1.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[23] == 1.0 );

    // Elements.
    TEST_ASSERT( mesh_block->elementIds()[0] == 12 );

    // Connectivity.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->connectivity()[i] == i );
    }

    // Permutation.
    for ( int i = 0; i < num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->permutation()[i] == i );
    }
}

//---------------------------------------------------------------------------//
// Pyramid mesh.
TEUCHOS_UNIT_TEST( MeshBlock, pyramid_block_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = buildPyramidBlock();

    // Mesh parameters.
    int vertex_dim = 3;
    int num_vertices = 5;
    int element_topo = DataTransferKit::DTK_PYRAMID;

    // Basic block info.
    TEST_ASSERT( (int) mesh_block->dimension() == vertex_dim );
    TEST_ASSERT( (int) mesh_block->verticesPerElement() == num_vertices );
    TEST_ASSERT( (int) mesh_block->elementTopology() == element_topo );

    // Vertices.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->vertexIds()[i] == i );
    }

    // Coords.
    // x
    TEST_ASSERT( mesh_block->vertexCoordinates()[0] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[1] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[2] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[3] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[4] == 0.5 );

    // y
    TEST_ASSERT( mesh_block->vertexCoordinates()[5]  == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[6]  == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[7] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[8] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[9] == 0.5 ); 

    // z
    TEST_ASSERT( mesh_block->vertexCoordinates()[10] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[11] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[12] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[13] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[14] == 1.0 );

    // Elements.
    TEST_ASSERT( mesh_block->elementIds()[0] == 12 );

    // Connectivity.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->connectivity()[i] == i );
    }

    // Permutation.
    for ( int i = 0; i < num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->permutation()[i] == i );
    }
}

//---------------------------------------------------------------------------//
// Wedge mesh.
TEUCHOS_UNIT_TEST( MeshBlock, wedge_block_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = buildWedgeBlock();

    // Mesh parameters.
    int vertex_dim = 3;
    int num_vertices = 6;
    int element_topo = DataTransferKit::DTK_WEDGE;

    // Basic block info.
    TEST_ASSERT( (int) mesh_block->dimension() == vertex_dim );
    TEST_ASSERT( (int) mesh_block->verticesPerElement() == num_vertices );
    TEST_ASSERT( (int) mesh_block->elementTopology() == element_topo );

    // Vertices.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->vertexIds()[i] == i );
    }

    // Coords.
    // x
    TEST_ASSERT( mesh_block->vertexCoordinates()[0] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[1] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[2] == 0.5 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[3] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[4] == 1.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[5] == 0.5 );

    // y
    TEST_ASSERT( mesh_block->vertexCoordinates()[6]  == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[7]  == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[8]  == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[9]  == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[10] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[11] == 1.0 ); 

    // z
    TEST_ASSERT( mesh_block->vertexCoordinates()[12] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[13] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[14] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[15] == 1.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[16] == 1.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[17] == 1.0 );

    // Elements.
    TEST_ASSERT( mesh_block->elementIds()[0] == 12 );

    // Connectivity.
    for ( unsigned i = 0; i < (unsigned) num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->connectivity()[i] == i );
    }

    // Permutation.
    for ( int i = 0; i < num_vertices; ++i )
    {
	TEST_ASSERT( mesh_block->permutation()[i] == i );
    }
}

//---------------------------------------------------------------------------//
// end tstMeshBlock.cpp
//---------------------------------------------------------------------------//
