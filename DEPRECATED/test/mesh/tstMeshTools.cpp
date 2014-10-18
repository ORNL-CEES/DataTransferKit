//---------------------------------------------------------------------------//
/*!
 * \file tstMeshTools.cpp
 * \author Stuart R. Slattery
 * \brief Mesh tools unit tests.
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

#include <DTK_MeshTools.hpp>
#include <DTK_MeshBlock.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_BoundingBox.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ScalarTraits.hpp>
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
TEUCHOS_UNIT_TEST( MeshBlock, line_tools_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildLineBlock();

    // Bounding Boxes.
    DataTransferKit::BoundingBox local_box = 
	DataTransferKit::MeshTools::localBoundingBox( mesh_block );
    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( local_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( local_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );

    DataTransferKit::BoundingBox global_box = 
	DataTransferKit::MeshTools::globalBoundingBox( mesh_block, 
						       getDefaultComm<int>() );
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
TEUCHOS_UNIT_TEST( MeshBlock, tri_tools_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildTriBlock();

    // Bounding Boxes.
    DataTransferKit::BoundingBox local_box = 
	DataTransferKit::MeshTools::localBoundingBox( mesh_block );
    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );

    DataTransferKit::BoundingBox global_box = 
	DataTransferKit::MeshTools::globalBoundingBox( mesh_block, 
						       getDefaultComm<int>() );
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
TEUCHOS_UNIT_TEST( MeshBlock, quad_tools_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildQuadBlock();

    // Bounding boxes.
    DataTransferKit::BoundingBox local_box = 
	DataTransferKit::MeshTools::localBoundingBox( mesh_block );
    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == -Teuchos::ScalarTraits<double>::rmax() );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == Teuchos::ScalarTraits<double>::rmax() );

    DataTransferKit::BoundingBox global_box = 
	DataTransferKit::MeshTools::globalBoundingBox( mesh_block, 
						       getDefaultComm<int>() );
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
TEUCHOS_UNIT_TEST( MeshBlock, tet_tools_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildTetBlock();

    // Bounding Boxes.
    DataTransferKit::BoundingBox local_box = 
	DataTransferKit::MeshTools::localBoundingBox( mesh_block );
    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == 0.0 );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == 1.0 );

    DataTransferKit::BoundingBox global_box = 
	DataTransferKit::MeshTools::globalBoundingBox( mesh_block, 
						       getDefaultComm<int>() );
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
TEUCHOS_UNIT_TEST( MeshBlock, hex_tools_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildHexBlock();

    // Bounding Boxes.
    DataTransferKit::BoundingBox local_box = 
	DataTransferKit::MeshTools::localBoundingBox( mesh_block );
    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == 0.0 );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == 1.0 );

    DataTransferKit::BoundingBox global_box = 
	DataTransferKit::MeshTools::globalBoundingBox( mesh_block, 
						       getDefaultComm<int>() );
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
TEUCHOS_UNIT_TEST( MeshBlock, pyramid_tools_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildPyramidBlock();

    // Bounding Boxes.
    DataTransferKit::BoundingBox local_box = 
	DataTransferKit::MeshTools::localBoundingBox( mesh_block );
    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == 0.0 );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == 1.0 );

    DataTransferKit::BoundingBox global_box = 
	DataTransferKit::MeshTools::globalBoundingBox( mesh_block, 
						       getDefaultComm<int>() );
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
TEUCHOS_UNIT_TEST( MeshBlock, wedge_tools_test )
{
    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildWedgeBlock();

    // Bounding Boxes.
    DataTransferKit::BoundingBox local_box = 
	DataTransferKit::MeshTools::localBoundingBox( mesh_block );
    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == 0.0 );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == 1.0 );

    DataTransferKit::BoundingBox global_box = 
	DataTransferKit::MeshTools::globalBoundingBox( mesh_block, 
						       getDefaultComm<int>() );
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
TEUCHOS_UNIT_TEST( MeshBlock, parallel_hex_tools_test )
{
    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Create a mesh block.
    Teuchos::RCP<DataTransferKit::MeshBlock> mesh_block = 
	DataTransferKit::UnitTest::MeshBuilder::buildParallelHexBlock();

    // Bounding Boxes.
    DataTransferKit::BoundingBox local_box = 
	DataTransferKit::MeshTools::localBoundingBox( mesh_block );
    Teuchos::Tuple<double,6> local_bounds = local_box.getBounds();
    TEST_ASSERT( local_bounds[0] == 0.0 );
    TEST_ASSERT( local_bounds[1] == 0.0 );
    TEST_ASSERT( local_bounds[2] == my_rank );
    TEST_ASSERT( local_bounds[3] == 1.0 );
    TEST_ASSERT( local_bounds[4] == 1.0 );
    TEST_ASSERT( local_bounds[5] == my_rank+1 );

    DataTransferKit::BoundingBox global_box = 
	DataTransferKit::MeshTools::globalBoundingBox( mesh_block, 
						       getDefaultComm<int>() );
    Teuchos::Tuple<double,6> global_bounds = global_box.getBounds();
    TEST_ASSERT( global_bounds[0] == 0.0 );
    TEST_ASSERT( global_bounds[1] == 0.0 );
    TEST_ASSERT( global_bounds[2] == 0.0 );
    TEST_ASSERT( global_bounds[3] == 1.0 );
    TEST_ASSERT( global_bounds[4] == 1.0 );
    TEST_ASSERT( global_bounds[5] == my_size );
}

//---------------------------------------------------------------------------//
// end tstMeshTools.cpp
//---------------------------------------------------------------------------//

