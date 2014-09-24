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
#include <DTK_MeshManager.hpp>
#include <DTK_MeshContainer.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>

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
TEUCHOS_UNIT_TEST( MeshContainer, line_element_tree_test )
{
    using namespace DataTransferKit;

    // Create a mesh container.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = 
	DataTransferKit::UnitTest::MeshBuilder::buildLineBlock();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager(mesh_blocks, getDefaultComm<int>(), 1) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(1);
    MeshId ordinal = 0;
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = 
	DataTransferKit::UnitTest::MeshBuilder::buildTriBlock();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager(mesh_blocks, getDefaultComm<int>(), 2) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(2);
    MeshId ordinal = 0;
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = 
	DataTransferKit::UnitTest::MeshBuilder::buildQuadBlock();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager(mesh_blocks, getDefaultComm<int>(), 2) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(2);
    MeshId ordinal = 0;
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = 
	DataTransferKit::UnitTest::MeshBuilder::buildTetBlock();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager(mesh_blocks, getDefaultComm<int>(), 3) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(3);
    MeshId ordinal = 0;
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = 
	DataTransferKit::UnitTest::MeshBuilder::buildHexBlock();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager(mesh_blocks, getDefaultComm<int>(), 3) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(3);
    MeshId ordinal = 0;
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = 
	DataTransferKit::UnitTest::MeshBuilder::buildPyramidBlock();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager(mesh_blocks, getDefaultComm<int>(), 3) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(3);
    MeshId ordinal = 0;

    point[0] = 0.25;
    point[1] = 0.25;
    point[2] = 0.25;
    TEST_ASSERT( element_tree.findPoint( point, ordinal ) );
    TEST_ASSERT( ordinal == 12 );

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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = 
	DataTransferKit::UnitTest::MeshBuilder::buildWedgeBlock();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager(mesh_blocks, getDefaultComm<int>(), 3) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(3);
    MeshId ordinal = 0;
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = 
	DataTransferKit::UnitTest::MeshBuilder::buildParallelHexBlock();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager(mesh_blocks, getDefaultComm<int>(), 3) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(3);

    MeshId ordinal = 0;
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 2 );
    mesh_blocks[0] = 
	DataTransferKit::UnitTest::MeshBuilder::buildTriBlock();
    mesh_blocks[1] = 
	DataTransferKit::UnitTest::MeshBuilder::buildShiftedQuadBlock();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager(mesh_blocks, getDefaultComm<int>(), 2) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree element_tree( mesh_manager );

    // Search the tree for some random points.
    int num_points = 1000;
    Teuchos::Array<double> point(2);
    MeshId ordinal = 0;
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
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 3 );
    mesh_blocks[0] = 
	DataTransferKit::UnitTest::MeshBuilder::buildTetBlock();
    mesh_blocks[1] = 
	DataTransferKit::UnitTest::MeshBuilder::buildShiftedPyramidBlock();
    mesh_blocks[2] = 
	DataTransferKit::UnitTest::MeshBuilder::buildShiftedHexBlock();

    // Create a mesh manager.
    Teuchos::RCP<MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager(mesh_blocks, getDefaultComm<int>(), 3) );
    mesh_manager->buildIndexing();

    // Create an element tree.
    ElementTree element_tree( mesh_manager );

    // Search the tree for some random points.
    double tol = 1.0e-8;
    int num_points = 1000;
    Teuchos::Array<double> point(3);
    MeshId ordinal = 0;
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
