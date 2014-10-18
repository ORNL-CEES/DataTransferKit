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

#include "DTK_TestMeshBuilder.hpp"

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Line mesh.
TEUCHOS_UNIT_TEST( MeshBlock, line_block_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildLineBlock();

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
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildTriBlock();

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
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildQuadBlock();

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
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildTetBlock();

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
    TEST_ASSERT( mesh_block->vertexCoordinates()[2] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[3] == 0.0 );

    // y
    TEST_ASSERT( mesh_block->vertexCoordinates()[4] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[5] == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[6] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[7] == 0.0 ); 

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
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildHexBlock();

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
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildPyramidBlock();

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
    TEST_ASSERT( mesh_block->vertexCoordinates()[4] == 0.0 );

    // y
    TEST_ASSERT( mesh_block->vertexCoordinates()[5]  == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[6]  == 0.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[7] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[8] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[9] == 0.0 ); 

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
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildWedgeBlock();

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
    TEST_ASSERT( mesh_block->vertexCoordinates()[2] == 1.0 ); 
    TEST_ASSERT( mesh_block->vertexCoordinates()[3] == 0.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[4] == 1.0 );
    TEST_ASSERT( mesh_block->vertexCoordinates()[5] == 1.0 );

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
