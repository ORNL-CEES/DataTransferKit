//---------------------------------------------------------------------------//
/*!
 * \file tstMeshManager.cpp
 * \author Stuart R. Slattery
 * \brief Mesh manager unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_MeshManager.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshBlock.hpp>
#include <DTK_MeshTools.hpp>
#include <DTK_BoundingBox.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>

#include "DTK_TestMeshBuilder.hpp"

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
// Tests
//---------------------------------------------------------------------------//
// Line mesh.
TEUCHOS_UNIT_TEST( MeshManager, line_manager_test )
{
    // Create a mesh container.
    typedef DataTransferKit::MeshBlock MeshType;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = DataTransferKit::UnitTest::MeshBuilder::buildLineBlock();

    // Create a mesh manager.
    DataTransferKit::MeshManager mesh_manager( mesh_blocks, getDefaultComm<int>(), 1 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 1 );

    // Mesh parameters.
    int num_vertices = 2;

    // Check the mesh data.
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> > blocks =
	mesh_manager.meshBlocks();
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> >::const_iterator block_iterator;
    for ( block_iterator = blocks.begin();
	  block_iterator != blocks.end();
	  ++block_iterator )
    {
	// Basic block info.
	TEST_ASSERT( (*block_iterator)->elementIds().size() == 1 );
	TEST_ASSERT( (*block_iterator)->vertexIds().size() == num_vertices );

	// Vertices.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> vertices_view =
	    (*block_iterator)->vertexIds();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    (*block_iterator)->vertexCoordinates();
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 

	// Elements.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> elements_view =
	    (*block_iterator)->elementIds();
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> connectivity_view =
	    (*block_iterator)->connectivity();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    (*block_iterator)->permutation();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding Boxes.
    DataTransferKit::BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );
}

//---------------------------------------------------------------------------//
// Tri mesh.
TEUCHOS_UNIT_TEST( MeshBlock, tri_manager_test )
{
    // Create a mesh block.
    typedef DataTransferKit::MeshBlock MeshType;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = DataTransferKit::UnitTest::MeshBuilder::buildTriBlock();

    // Create a mesh manager.
    DataTransferKit::MeshManager mesh_manager( mesh_blocks, getDefaultComm<int>(), 2 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 2 );

    // Mesh parameters.
    int num_vertices = 3;

    // Check the mesh data.
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> > blocks =
	mesh_manager.meshBlocks();
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> >::const_iterator block_iterator;
    for ( block_iterator = blocks.begin();
	  block_iterator != blocks.end();
	  ++block_iterator )
    {
	TEST_EQUALITY( (*block_iterator)->vertexIds().size(), num_vertices );
	TEST_EQUALITY( (*block_iterator)->elementIds().size(), 1 );

	// Vertices.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> vertices_view = 
	    (*block_iterator)->vertexIds();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    (*block_iterator)->vertexCoordinates();
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 1.0 ); 

	// y
	TEST_ASSERT( coords_view[3] == 0.0 ); 
	TEST_ASSERT( coords_view[4] == 0.0 ); 
	TEST_ASSERT( coords_view[5] == 1.0 ); 

	// Elements.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> elements_view =
	    (*block_iterator)->elementIds();
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> connectivity_view =
	    (*block_iterator)->connectivity();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    (*block_iterator)->permutation();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    DataTransferKit::BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );
}

//---------------------------------------------------------------------------//
// Quad mesh.
TEUCHOS_UNIT_TEST( MeshBlock, quad_manager_test )
{
    // Create a mesh block.
    typedef DataTransferKit::MeshBlock MeshType;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = DataTransferKit::UnitTest::MeshBuilder::buildQuadBlock();

    // Create a mesh manager.
    DataTransferKit::MeshManager mesh_manager( mesh_blocks, getDefaultComm<int>(), 2 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 2 );

    // Mesh parameters.
    int num_vertices = 4;

    // Check the mesh data.
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> > blocks =
	mesh_manager.meshBlocks();
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> >::const_iterator block_iterator;
    for ( block_iterator = blocks.begin();
	  block_iterator != blocks.end();
	  ++block_iterator )
    {
	// Vertices.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> vertices_view = 
	    (*block_iterator)->vertexIds();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    (*block_iterator)->vertexCoordinates();
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 1.0 ); 
	TEST_ASSERT( coords_view[3] == 0.0 );

	// y
	TEST_ASSERT( coords_view[4] == 0.0 ); 
	TEST_ASSERT( coords_view[5] == 0.0 ); 
	TEST_ASSERT( coords_view[6] == 1.0 ); 
	TEST_ASSERT( coords_view[7] == 1.0 ); 

	// Elements.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> elements_view =
	    (*block_iterator)->elementIds();
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> connectivity_view =
	    (*block_iterator)->connectivity();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    (*block_iterator)->permutation();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding boxes.
    DataTransferKit::BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );
}

//---------------------------------------------------------------------------//
// Tet mesh.
TEUCHOS_UNIT_TEST( MeshBlock, tet_manager_test )
{
    // Create a mesh block.
    typedef DataTransferKit::MeshBlock MeshType;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = DataTransferKit::UnitTest::MeshBuilder::buildTetBlock();

    // Create a mesh manager.
    DataTransferKit::MeshManager mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Mesh parameters.
    int num_vertices = 4;

    // Check the mesh data.
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> > blocks =
	mesh_manager.meshBlocks();
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> >::const_iterator block_iterator;
    for ( block_iterator = blocks.begin();
	  block_iterator != blocks.end();
	  ++block_iterator )
    {
	// Vertices.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> vertices_view = 
	    (*block_iterator)->vertexIds();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    (*block_iterator)->vertexCoordinates();
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 0.0 ); 
	TEST_ASSERT( coords_view[3] == 0.0 );

	// y
	TEST_ASSERT( coords_view[4] == 0.0 ); 
	TEST_ASSERT( coords_view[5] == 0.0 ); 
	TEST_ASSERT( coords_view[6] == 1.0 ); 
	TEST_ASSERT( coords_view[7] == 0.0 ); 

	// z
	TEST_ASSERT( coords_view[8]  == 0.0 );
	TEST_ASSERT( coords_view[9]  == 0.0 );
	TEST_ASSERT( coords_view[10] == 0.0 );
	TEST_ASSERT( coords_view[11] == 1.0 );

	// Elements.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> elements_view =
	    (*block_iterator)->elementIds();
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> connectivity_view =
	    (*block_iterator)->connectivity();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    (*block_iterator)->permutation();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding Boxes.
    DataTransferKit::BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// Hex mesh.
TEUCHOS_UNIT_TEST( MeshBlock, hex_manager_test )
{
    // Create a mesh block.
    typedef DataTransferKit::MeshBlock MeshType;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = DataTransferKit::UnitTest::MeshBuilder::buildHexBlock();

    // Create a mesh manager.
    DataTransferKit::MeshManager mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Mesh parameters.
    int num_vertices = 8;

    // Check the mesh data.
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> > blocks =
	mesh_manager.meshBlocks();
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> >::const_iterator block_iterator;
    for ( block_iterator = blocks.begin();
	  block_iterator != blocks.end();
	  ++block_iterator )
    {
	// Vertices.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> vertices_view = 
	    (*block_iterator)->vertexIds();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    (*block_iterator)->vertexCoordinates();
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 1.0 ); 
	TEST_ASSERT( coords_view[3] == 0.0 );
	TEST_ASSERT( coords_view[4] == 0.0 );
	TEST_ASSERT( coords_view[5] == 1.0 ); 
	TEST_ASSERT( coords_view[6] == 1.0 ); 
	TEST_ASSERT( coords_view[7] == 0.0 ); 

	// y
	TEST_ASSERT( coords_view[8]  == 0.0 ); 
	TEST_ASSERT( coords_view[9]  == 0.0 ); 
	TEST_ASSERT( coords_view[10] == 1.0 ); 
	TEST_ASSERT( coords_view[11] == 1.0 ); 
	TEST_ASSERT( coords_view[12] == 0.0 ); 
	TEST_ASSERT( coords_view[13] == 0.0 );
	TEST_ASSERT( coords_view[14] == 1.0 );
	TEST_ASSERT( coords_view[15] == 1.0 );

	// z
	TEST_ASSERT( coords_view[16] == 0.0 );
	TEST_ASSERT( coords_view[17] == 0.0 );
	TEST_ASSERT( coords_view[18] == 0.0 );
	TEST_ASSERT( coords_view[19] == 0.0 );
	TEST_ASSERT( coords_view[20] == 1.0 );
	TEST_ASSERT( coords_view[21] == 1.0 );
	TEST_ASSERT( coords_view[22] == 1.0 );
	TEST_ASSERT( coords_view[23] == 1.0 );

	// Elements.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> elements_view =
	    (*block_iterator)->elementIds();
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> connectivity_view =
	    (*block_iterator)->connectivity();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    (*block_iterator)->permutation();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding Boxes.
    DataTransferKit::BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// Pyramid mesh.
TEUCHOS_UNIT_TEST( MeshBlock, pyramid_manager_test )
{
    // Create a mesh block.
    typedef DataTransferKit::MeshBlock MeshType;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = DataTransferKit::UnitTest::MeshBuilder::buildPyramidBlock();

    // Create a mesh manager.
    DataTransferKit::MeshManager mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Mesh parameters.
    int num_vertices = 5;

    // Check the mesh data.
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> > blocks =
	mesh_manager.meshBlocks();
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> >::const_iterator block_iterator;
    for ( block_iterator = blocks.begin();
	  block_iterator != blocks.end();
	  ++block_iterator )
    {
	// Vertices.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> vertices_view = 
	    (*block_iterator)->vertexIds();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    (*block_iterator)->vertexCoordinates();
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 1.0 ); 
	TEST_ASSERT( coords_view[3] == 0.0 );
	TEST_ASSERT( coords_view[4] == 0.0 );

	// y
	TEST_ASSERT( coords_view[5]  == 0.0 ); 
	TEST_ASSERT( coords_view[6]  == 0.0 ); 
	TEST_ASSERT( coords_view[7] == 1.0 ); 
	TEST_ASSERT( coords_view[8] == 1.0 ); 
	TEST_ASSERT( coords_view[9] == 0.0 ); 

	// z
	TEST_ASSERT( coords_view[10] == 0.0 );
	TEST_ASSERT( coords_view[11] == 0.0 );
	TEST_ASSERT( coords_view[12] == 0.0 );
	TEST_ASSERT( coords_view[13] == 0.0 );
	TEST_ASSERT( coords_view[14] == 1.0 );

	// Elements.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> elements_view =
	    (*block_iterator)->elementIds();
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> connectivity_view =
	    (*block_iterator)->connectivity();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    (*block_iterator)->permutation();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding Boxes.
    DataTransferKit::BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// Wedge mesh.
TEUCHOS_UNIT_TEST( MeshBlock, wedge_manager_test )
{
    // Create a mesh block.
    typedef DataTransferKit::MeshBlock MeshType;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = DataTransferKit::UnitTest::MeshBuilder::buildWedgeBlock();

    // Create a mesh manager.
    DataTransferKit::MeshManager mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Mesh parameters.
    int num_vertices = 6;

    // Check the mesh data.
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> > blocks =
	mesh_manager.meshBlocks();
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> >::const_iterator block_iterator;
    for ( block_iterator = blocks.begin();
	  block_iterator != blocks.end();
	  ++block_iterator )
    {
	// Vertices.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> vertices_view = 
	    (*block_iterator)->vertexIds();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    (*block_iterator)->vertexCoordinates();
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 1.0 ); 
	TEST_ASSERT( coords_view[3] == 0.0 );
	TEST_ASSERT( coords_view[4] == 1.0 );
	TEST_ASSERT( coords_view[5] == 1.0 );

	// y
	TEST_ASSERT( coords_view[6]  == 0.0 ); 
	TEST_ASSERT( coords_view[7]  == 0.0 ); 
	TEST_ASSERT( coords_view[8]  == 1.0 ); 
	TEST_ASSERT( coords_view[9]  == 0.0 ); 
	TEST_ASSERT( coords_view[10] == 0.0 ); 
	TEST_ASSERT( coords_view[11] == 1.0 ); 

	// z
	TEST_ASSERT( coords_view[12] == 0.0 );
	TEST_ASSERT( coords_view[13] == 0.0 );
	TEST_ASSERT( coords_view[14] == 0.0 );
	TEST_ASSERT( coords_view[15] == 1.0 );
	TEST_ASSERT( coords_view[16] == 1.0 );
	TEST_ASSERT( coords_view[17] == 1.0 );

	// Elements.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> elements_view =
	    (*block_iterator)->elementIds();
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> connectivity_view =
	    (*block_iterator)->connectivity();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    (*block_iterator)->permutation();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding Boxes.
    DataTransferKit::BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// Parallel hex mesh.
TEUCHOS_UNIT_TEST( MeshBlock, parallel_hex_manager_test )
{
    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Create a mesh block.
    typedef DataTransferKit::MeshBlock MeshType;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 1 );
    mesh_blocks[0] = DataTransferKit::UnitTest::MeshBuilder::buildParallelHexBlock();

    // Create a mesh manager.
    DataTransferKit::MeshManager mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 1 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Mesh parameters.
    int num_vertices = 8;

    // Check the mesh data.
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> > blocks =
	mesh_manager.meshBlocks();
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> >::const_iterator block_iterator;
    for ( block_iterator = blocks.begin();
	  block_iterator != blocks.end();
	  ++block_iterator )
    {
	// Vertices.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> vertices_view = 
	    (*block_iterator)->vertexIds();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Coords.
	Teuchos::ArrayRCP<const double> coords_view = 
	    (*block_iterator)->vertexCoordinates();
	// x
	TEST_ASSERT( coords_view[0] == 0.0 ); 
	TEST_ASSERT( coords_view[1] == 1.0 ); 
	TEST_ASSERT( coords_view[2] == 1.0 ); 
	TEST_ASSERT( coords_view[3] == 0.0 );
	TEST_ASSERT( coords_view[4] == 0.0 );
	TEST_ASSERT( coords_view[5] == 1.0 ); 
	TEST_ASSERT( coords_view[6] == 1.0 ); 
	TEST_ASSERT( coords_view[7] == 0.0 ); 

	// y
	TEST_ASSERT( coords_view[8]  == 0.0 ); 
	TEST_ASSERT( coords_view[9]  == 0.0 ); 
	TEST_ASSERT( coords_view[10] == 1.0 ); 
	TEST_ASSERT( coords_view[11] == 1.0 ); 
	TEST_ASSERT( coords_view[12] == 0.0 ); 
	TEST_ASSERT( coords_view[13] == 0.0 );
	TEST_ASSERT( coords_view[14] == 1.0 );
	TEST_ASSERT( coords_view[15] == 1.0 );

	// z
	TEST_ASSERT( coords_view[16] == my_rank );
	TEST_ASSERT( coords_view[17] == my_rank );
	TEST_ASSERT( coords_view[18] == my_rank );
	TEST_ASSERT( coords_view[19] == my_rank );
	TEST_ASSERT( coords_view[20] == my_rank+1 );
	TEST_ASSERT( coords_view[21] == my_rank+1 );
	TEST_ASSERT( coords_view[22] == my_rank+1 );
	TEST_ASSERT( coords_view[23] == my_rank+1 );

	// Elements.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> elements_view =
	    (*block_iterator)->elementIds();
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> connectivity_view =
	    (*block_iterator)->connectivity();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    (*block_iterator)->permutation();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    // Bounding Boxes.
    DataTransferKit::BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == my_size );
}

//---------------------------------------------------------------------------//
// 2d hybrid test.
TEUCHOS_UNIT_TEST( MeshBlock, 2d_hybrid_manager_test )
{
    // Create a mesh block.
    typedef DataTransferKit::MeshBlock MeshType;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 2 );
    mesh_blocks[0] = DataTransferKit::UnitTest::MeshBuilder::buildTriBlock();
    mesh_blocks[1] = DataTransferKit::UnitTest::MeshBuilder::buildQuadBlock();

    // Create a mesh manager.
    DataTransferKit::MeshManager mesh_manager( mesh_blocks, getDefaultComm<int>(), 2 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 2 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );

    TEST_ASSERT( mesh_manager.dim() == 2 );

    // Check the mesh data.
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> > blocks =
	mesh_manager.meshBlocks();
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> >::const_iterator block_iterator;
    for ( block_iterator = blocks.begin();
	  block_iterator != blocks.end();
	  ++block_iterator )
    {
	TEST_ASSERT( (*block_iterator)->elementIds().size() == 1 );
	int num_vertices = (*block_iterator)->vertexIds().size();

	// Vertices.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> vertices_view = 
	    (*block_iterator)->vertexIds();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Elements.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> elements_view =
	    (*block_iterator)->elementIds();
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> connectivity_view =
	    (*block_iterator)->connectivity();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    (*block_iterator)->permutation();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    DataTransferKit::BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );
}

//---------------------------------------------------------------------------//
// 3d hybrid test.
TEUCHOS_UNIT_TEST( MeshBlock, 3d_hybrid_manager_test )
{
    // Create a mesh block.
    typedef DataTransferKit::MeshBlock MeshType;
    Teuchos::ArrayRCP<Teuchos::RCP<MeshType> > mesh_blocks( 3 );
    mesh_blocks[0] = DataTransferKit::UnitTest::MeshBuilder::buildTetBlock();
    mesh_blocks[1] = DataTransferKit::UnitTest::MeshBuilder::buildHexBlock();
    mesh_blocks[2] = DataTransferKit::UnitTest::MeshBuilder::buildPyramidBlock();

    // Create a mesh manager.
    DataTransferKit::MeshManager mesh_manager( mesh_blocks, getDefaultComm<int>(), 3 );
    TEST_ASSERT( mesh_manager.getNumBlocks() == 3 );
    TEST_ASSERT( mesh_manager.comm()->getRank() == getDefaultComm<int>()->getRank() );
    TEST_ASSERT( mesh_manager.comm()->getSize() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( mesh_manager.dim() == 3 );

    // Check the mesh data.
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> > blocks =
	mesh_manager.meshBlocks();
    Teuchos::ArrayRCP<Teuchos::RCP<DataTransferKit::MeshBlock> >::const_iterator block_iterator;
    for ( block_iterator = blocks.begin();
	  block_iterator != blocks.end();
	  ++block_iterator )
    {
	TEST_ASSERT( (*block_iterator)->elementIds().size() == 1 );
	int num_vertices = (*block_iterator)->vertexIds().size();

	// Vertices.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> vertices_view = 
	    (*block_iterator)->vertexIds();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) vertices_view[i] == i );
	}

	// Elements.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> elements_view =
	    (*block_iterator)->elementIds();
	TEST_ASSERT( elements_view[0] == 12 );

	// Connectivity.
	Teuchos::ArrayRCP<const DataTransferKit::MeshId> connectivity_view =
	    (*block_iterator)->connectivity();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) connectivity_view[i] == i );
	}

	// Permutation.
	Teuchos::ArrayRCP<const int> permutation_view =
	    (*block_iterator)->permutation();
	for ( int i = 0; i < num_vertices; ++i )
	{
	    TEST_ASSERT( (int) permutation_view[i] == i );
	}
    }

    DataTransferKit::BoundingBox global_box = mesh_manager.globalBoundingBox();
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == 1.0 );
}

//---------------------------------------------------------------------------//
// end tstMeshManager.cpp
//---------------------------------------------------------------------------//

