//---------------------------------------------------------------------------//
/*!
 * \file tstRCB.cpp
 * \author Stuart R. Slattery
 * \brief Unit tests for recursive coordinate bisectioning.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_RCB.hpp>
#include <DTK_BoundingBox.hpp>
#include <DTK_MeshTypes.hpp>
#include <DTK_MeshManager.hpp>
#include <DTK_MeshBlock.hpp>
#include <DTK_MeshContainer.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
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
    return Teuchos::rcp( new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Global test variables.
//---------------------------------------------------------------------------//
int rand_size = 1000;

//---------------------------------------------------------------------------//
// Mesh container creation functions.
//---------------------------------------------------------------------------//
// 1d mesh
Teuchos::RCP<DataTransferKit::MeshBlock> build1dContainer()
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();

    // Make some random vertices.
    int vertex_dim = 1;
    int num_rand = vertex_dim*rand_size;
    std::srand( 1 );
    Teuchos::ArrayRCP<MeshId> vertex_handles( rand_size );
    Teuchos::ArrayRCP<double> coords( num_rand );
    for ( int i = 0; i < rand_size; ++i )
    {
	vertex_handles[i] = i*rand_size + my_rank;
	coords[i] = (double) std::rand() / RAND_MAX;
    }

    // Empty element vectors. We only need vertices for these tests.
    Teuchos::ArrayRCP<MeshId> element_handles;
    Teuchos::ArrayRCP<MeshId> element_connectivity;

    Teuchos::ArrayRCP<int> permutation_list(2);
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i; 
    }

    return Teuchos::rcp( 
	new MeshContainer( vertex_dim, vertex_handles, coords, 
			   DTK_LINE_SEGMENT, 2,
			   element_handles, element_connectivity,
			   permutation_list ) );
}

//---------------------------------------------------------------------------//
// 2d mesh
Teuchos::RCP<DataTransferKit::MeshBlock> build2dContainer()
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();

    // Make some random vertices.
    int vertex_dim = 2;
    int num_rand = vertex_dim*rand_size;
    std::srand( 1 );
    Teuchos::ArrayRCP<MeshId> vertex_handles( rand_size );
    Teuchos::ArrayRCP<double> coords( num_rand );
    for ( int i = 0; i < rand_size; ++i )
    {
	vertex_handles[i] = i*rand_size + my_rank;
	coords[ i ] = (double) std::rand() / RAND_MAX;
	coords[ rand_size + i ] = (double) std::rand() / RAND_MAX;
    }

    // Empty element vectors. We only need vertices for these tests.
    Teuchos::ArrayRCP<MeshId> element_handles;
    Teuchos::ArrayRCP<MeshId> element_connectivity;
    Teuchos::ArrayRCP<int> permutation_list(3);
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i; 
    }

    return Teuchos::rcp( 
	new MeshContainer( vertex_dim, vertex_handles, coords, 
			   DTK_TRIANGLE, 3,
			   element_handles, element_connectivity,
			   permutation_list ) );
}

//---------------------------------------------------------------------------//
// 3d mesh
Teuchos::RCP<DataTransferKit::MeshBlock> build3dContainer()
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();

    // Make some random vertices.
    int vertex_dim = 3;
    int num_rand = vertex_dim*rand_size;
    std::srand( 1 );
    Teuchos::ArrayRCP<MeshId> vertex_handles( rand_size );
    Teuchos::ArrayRCP<double> coords( num_rand );
    for ( int i = 0; i < rand_size; ++i )
    {
	vertex_handles[i] = i*rand_size + my_rank;
	coords[ i ] = (double) std::rand() / RAND_MAX;
	coords[ rand_size + i ] = (double) std::rand() / RAND_MAX;
	coords[ 2*rand_size + i ] = (double) std::rand() / RAND_MAX;
    }

    // Empty element vectors. We only need vertices for these tests.
    Teuchos::ArrayRCP<MeshId> element_handles;
    Teuchos::ArrayRCP<MeshId> element_connectivity;
    Teuchos::ArrayRCP<int> permutation_list(8);
    for ( int i = 0; i < permutation_list.size(); ++i )
    {
	permutation_list[i] = i; 
    }

    return Teuchos::rcp( 
	new MeshContainer( vertex_dim, vertex_handles, coords, 
			   DTK_HEXAHEDRON, 8,
			   element_handles, element_connectivity,
			   permutation_list ) );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// 1d mesh
TEUCHOS_UNIT_TEST( RCB, 1d_rcb_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // create a mesh.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = build1dContainer();

    // All of the vertices will be partitioned.
    int mesh_dim = 1;
    int num_vertices = mesh_blocks[0]->vertexIds().size();
    int num_coords = mesh_dim * num_vertices;
    Teuchos::Array<short int> active_vertices( num_vertices, 1 );

    // Create a mesh manager.
    Teuchos::RCP< MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager( 
	    mesh_blocks, getDefaultComm<int>(), mesh_dim ) );
    mesh_manager->setActiveVertices( active_vertices, 0 );

    // Partition the mesh with RCB.
    typedef RCB::zoltan_id_type zoltan_id_type;
    RCB rcb( getDefaultComm<int>(), mesh_manager, 1  );
    rcb.partition( mesh_manager->globalBoundingBox() );
    
    // Get the random numbers that were used to compute the vertex coordinates.
    std::srand( 1 );
    Teuchos::Array<double> random_numbers;
    for ( int i = 0; i < num_coords; ++i )
    {
	random_numbers.push_back( (double) std::rand() / RAND_MAX );
    }

    // Check that these are in fact the random numbers used for the vertices.
    Teuchos::ArrayRCP<const double> vertex_coords =
	mesh_blocks[0]->vertexCoordinates();
    for ( int i = 0; i < num_vertices; ++i )
    {
	TEST_ASSERT( vertex_coords[ i ] == 
		     random_numbers[mesh_dim*i] );
    }

    // Check import parameters.
    int num_import = rcb.getNumImport();
    TEST_ASSERT( num_import == rand_size*(my_size-1) / my_size );

    Teuchos::ArrayView<zoltan_id_type> import_global_ids = 
	rcb.getImportGlobalIds();
    TEST_ASSERT( import_global_ids.size() == num_import );

    Teuchos::ArrayView<zoltan_id_type> import_local_ids = 
	rcb.getImportLocalIds();
    TEST_ASSERT( import_local_ids.size() == num_import );

    Teuchos::ArrayView<int> import_procs = rcb.getImportProcs();
    TEST_ASSERT( import_procs.size() == num_import );

    Teuchos::ArrayView<int> import_parts = rcb.getImportParts();
    TEST_ASSERT( import_parts.size() == num_import );

    for ( int i = 0; i < num_import; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( import_procs[i] != my_rank &&
		     import_procs[i] >= 0 &&
		     import_procs[i] < my_size );

	TEST_ASSERT( import_parts[i] == my_rank );
    }

    // Check export parameters.
    int num_export = rcb.getNumExport();
    TEST_ASSERT( num_export == rand_size*(my_size-1) / my_size );

    Teuchos::ArrayView<zoltan_id_type> export_global_ids = 
	rcb.getExportGlobalIds();
    TEST_ASSERT( export_global_ids.size() == num_export );

    Teuchos::ArrayView<zoltan_id_type> export_local_ids = 
	rcb.getExportLocalIds();
    TEST_ASSERT( export_local_ids.size() == num_export );

    Teuchos::ArrayView<int> export_procs = rcb.getExportProcs();
    TEST_ASSERT( export_procs.size() == num_export );
    
    Teuchos::ArrayView<int> export_parts = rcb.getExportParts();
    TEST_ASSERT( export_parts.size() == num_export );

    for ( int i = 0; i < num_export; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( export_procs[i] == export_parts[i] );

	TEST_ASSERT( export_procs[i] != my_rank &&
		     export_procs[i] >= 0 &&
		     export_procs[i] < my_size );

	TEST_ASSERT( export_parts[i] != my_rank &&
		     export_parts[i] >= 0 &&
		     export_parts[i] < my_size );
    }

    // Check the destination proc point search.
    Teuchos::Array<double> point_0(1);
    point_0[0] = 2.0;
    TEST_ASSERT( rcb.getPointDestinationProc( point_0() ) == my_size-1 );

    Teuchos::Array<double> point_1(1);
    point_1[0] = -2.0;
    TEST_ASSERT( rcb.getPointDestinationProc( point_1() ) == 0 );

    Teuchos::Array<double> point_2(1);
    point_2[0] = 0.2;
    TEST_ASSERT( rcb.getPointDestinationProc( point_2() ) == 0 );

    Teuchos::Array<double> point_3(1);
    point_3[0] = 0.8;
    TEST_ASSERT( rcb.getPointDestinationProc( point_3() ) == my_size-1 );
}

//---------------------------------------------------------------------------//
// 2d mesh
TEUCHOS_UNIT_TEST( RCB, 2d_rcb_test )
{
    using namespace DataTransferKit;

    // create a mesh.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = build2dContainer();

    // All of the vertices will be partitioned.
    int mesh_dim = 2;
    int num_vertices = mesh_blocks[0]->vertexIds().size();
    int num_coords = mesh_dim * num_vertices;
    Teuchos::Array<short int> active_vertices( num_vertices, 1 );

    // Create a mesh manager.
    Teuchos::RCP<MeshManager> mesh_manager = Teuchos::rcp( 
	new MeshManager( 
	    mesh_blocks, getDefaultComm<int>(), mesh_dim ) );
    mesh_manager->setActiveVertices( active_vertices, 0 );

    // Partition the mesh with RCB.
    typedef RCB::zoltan_id_type zoltan_id_type;
    RCB rcb( getDefaultComm<int>(), mesh_manager, 2 );
    rcb.partition( mesh_manager->globalBoundingBox() );
    
    // Get the random numbers that were used to compute the vertex coordinates.
    std::srand( 1 );
    Teuchos::Array<double> random_numbers;
    for ( int i = 0; i < num_coords; ++i )
    {
	random_numbers.push_back( (double) std::rand() / RAND_MAX );
    }

    // Check that these are in fact the random numbers used for the vertices.
    Teuchos::ArrayRCP<const double> vertex_coords =
	mesh_blocks[0]->vertexCoordinates();
    for ( int i = 0; i < num_vertices; ++i )
    {
	TEST_ASSERT( vertex_coords[ i ] == 
		     random_numbers[mesh_dim*i] );
	TEST_ASSERT( vertex_coords[ num_vertices + i ] == 
		     random_numbers[mesh_dim*i+1] );
    }

    // Get MPI parameters.
    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Check import parameters.
    int num_import = rcb.getNumImport();
    TEST_ASSERT( num_import == rand_size*(my_size-1) / my_size );

    Teuchos::ArrayView<zoltan_id_type> import_global_ids = 
	rcb.getImportGlobalIds();
    TEST_ASSERT( import_global_ids.size() == num_import );

    Teuchos::ArrayView<zoltan_id_type> import_local_ids = 
	rcb.getImportLocalIds();
    TEST_ASSERT( import_local_ids.size() == num_import );

    Teuchos::ArrayView<int> import_procs = rcb.getImportProcs();
    TEST_ASSERT( import_procs.size() == num_import );

    Teuchos::ArrayView<int> import_parts = rcb.getImportParts();
    TEST_ASSERT( import_parts.size() == num_import );

    for ( int i = 0; i < num_import; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( import_procs[i] != my_rank &&
		     import_procs[i] >= 0 &&
		     import_procs[i] < my_size );

	TEST_ASSERT( import_parts[i] == my_rank );
    }

    // Check export parameters.
    int num_export = rcb.getNumExport();
    TEST_ASSERT( num_export == rand_size*(my_size-1) / my_size );

    Teuchos::ArrayView<zoltan_id_type> export_global_ids = 
	rcb.getExportGlobalIds();
    TEST_ASSERT( export_global_ids.size() == num_export );

    Teuchos::ArrayView<zoltan_id_type> export_local_ids = 
	rcb.getExportLocalIds();
    TEST_ASSERT( export_local_ids.size() == num_export );

    Teuchos::ArrayView<int> export_procs = rcb.getExportProcs();
    TEST_ASSERT( export_procs.size() == num_export );
    
    Teuchos::ArrayView<int> export_parts = rcb.getExportParts();
    TEST_ASSERT( export_parts.size() == num_export );

    for ( int i = 0; i < num_export; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( export_procs[i] == export_parts[i] );

	TEST_ASSERT( export_procs[i] != my_rank &&
		     export_procs[i] >= 0 &&
		     export_procs[i] < my_size );

	TEST_ASSERT( export_parts[i] != my_rank &&
		     export_parts[i] >= 0 &&
		     export_parts[i] < my_size );
    }

    // Check the destination proc point search.
    Teuchos::Array<double> point_0(2);
    point_0[0] = 2.0;
    point_0[1] = 2.0;
    TEST_ASSERT( rcb.getPointDestinationProc( point_0() ) == my_size-1 );

    Teuchos::Array<double> point_1(2);
    point_1[0] = -2.0;
    point_1[1] = -2.0;
    TEST_ASSERT( rcb.getPointDestinationProc( point_1() ) == 0 );

    Teuchos::Array<double> point_2(2);
    point_2[0] = 0.2;
    point_2[1] = 0.2;
    TEST_ASSERT( rcb.getPointDestinationProc( point_2() ) == 0 );

    Teuchos::Array<double> point_3(2);
    point_3[0] = 0.8;
    point_3[1] = 0.8;
    TEST_ASSERT( rcb.getPointDestinationProc( point_3() ) == my_size-1 );
}

//---------------------------------------------------------------------------//
// 3d mesh
TEUCHOS_UNIT_TEST( RCB, 3d_rcb_test )
{
    using namespace DataTransferKit;

    // create a mesh.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = build3dContainer();

    // All of the vertices will be partitioned.
    int mesh_dim = 3;
    int num_vertices = mesh_blocks[0]->vertexIds().size();
    int num_coords = mesh_dim * num_vertices;
    Teuchos::Array<short int> active_vertices( num_vertices, 1 );

    // Create a mesh manager.
    Teuchos::RCP< MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager( 
	    mesh_blocks, getDefaultComm<int>(), mesh_dim ) );
    mesh_manager->setActiveVertices( active_vertices, 0 );

    // Partition the mesh with RCB.
    typedef RCB::zoltan_id_type zoltan_id_type;
    RCB rcb( getDefaultComm<int>(), mesh_manager, 3 );
    rcb.partition( mesh_manager->globalBoundingBox() );
    
    // Get the random numbers that were used to compute the vertex coordinates.
    std::srand( 1 );
    Teuchos::Array<double> random_numbers;
    for ( int i = 0; i < num_coords; ++i )
    {
	random_numbers.push_back( (double) std::rand() / RAND_MAX );
    }

    // Check that these are in fact the random numbers used for the vertices.
    Teuchos::ArrayRCP<const double> vertex_coords = 
	mesh_blocks[0]->vertexCoordinates();
    for ( int i = 0; i < num_vertices; ++i )
    {
	TEST_ASSERT( vertex_coords[ i ] == 
		     random_numbers[mesh_dim*i] );
	TEST_ASSERT( vertex_coords[ num_vertices + i ] == 
		     random_numbers[mesh_dim*i+1] );
	TEST_ASSERT( vertex_coords[ 2*num_vertices + i ] == 
		     random_numbers[mesh_dim*i+2] );
    }

    // Get MPI parameters.
    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Check import parameters.
    int num_import = rcb.getNumImport();
    TEST_ASSERT( num_import == rand_size*(my_size-1) / my_size );

    Teuchos::ArrayView<zoltan_id_type> import_global_ids = 
	rcb.getImportGlobalIds();
    TEST_ASSERT( import_global_ids.size() == num_import );

    Teuchos::ArrayView<zoltan_id_type> import_local_ids = 
	rcb.getImportLocalIds();
    TEST_ASSERT( import_local_ids.size() == num_import );

    Teuchos::ArrayView<int> import_procs = rcb.getImportProcs();
    TEST_ASSERT( import_procs.size() == num_import );

    Teuchos::ArrayView<int> import_parts = rcb.getImportParts();
    TEST_ASSERT( import_parts.size() == num_import );

    for ( int i = 0; i < num_import; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( import_procs[i] != my_rank &&
		     import_procs[i] >= 0 &&
		     import_procs[i] < my_size );

	TEST_ASSERT( import_parts[i] == my_rank );
    }

    // Check export parameters.
    int num_export = rcb.getNumExport();
    TEST_ASSERT( num_export == rand_size*(my_size-1) / my_size );

    Teuchos::ArrayView<zoltan_id_type> export_global_ids = 
	rcb.getExportGlobalIds();
    TEST_ASSERT( export_global_ids.size() == num_export );

    Teuchos::ArrayView<zoltan_id_type> export_local_ids = 
	rcb.getExportLocalIds();
    TEST_ASSERT( export_local_ids.size() == num_export );

    Teuchos::ArrayView<int> export_procs = rcb.getExportProcs();
    TEST_ASSERT( export_procs.size() == num_export );
    
    Teuchos::ArrayView<int> export_parts = rcb.getExportParts();
    TEST_ASSERT( export_parts.size() == num_export );

    for ( int i = 0; i < num_export; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( export_procs[i] == export_parts[i] );

	TEST_ASSERT( export_procs[i] != my_rank &&
		     export_procs[i] >= 0 &&
		     export_procs[i] < my_size );

	TEST_ASSERT( export_parts[i] != my_rank &&
		     export_parts[i] >= 0 &&
		     export_parts[i] < my_size );
    }

    // Check the destination proc point search.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 2.0;
    point_0[1] = 2.0;
    point_0[2] = 2.0;
    TEST_ASSERT( rcb.getPointDestinationProc( point_0() ) == my_size-1 );

    Teuchos::Array<double> point_1(3);
    point_1[0] = -2.0;
    point_1[1] = -2.0;
    point_1[2] = -2.0;
    TEST_ASSERT( rcb.getPointDestinationProc( point_1() ) == 0 );

    Teuchos::Array<double> point_2(3);
    point_2[0] = 0.2;
    point_2[1] = 0.2;
    point_2[2] = 0.2;
    TEST_ASSERT( rcb.getPointDestinationProc( point_2() ) == 0 );

    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.8;
    point_3[1] = 0.8;
    point_3[2] = 0.8;
    TEST_ASSERT( rcb.getPointDestinationProc( point_3() ) == my_size-1 );
}

//---------------------------------------------------------------------------//
// partial 1d mesh
TEUCHOS_UNIT_TEST( RCB, partial_1d_rcb_test )
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // create a mesh.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = build1dContainer();

    // Only some of the vertices will be partitioned.
    int mesh_dim = 1;
    int num_vertices = mesh_blocks[0]->vertexIds().size();
    int num_coords = mesh_dim * num_vertices;
    Teuchos::Array<short int> active_vertices( num_vertices, 1 );
    active_vertices[0] = 0;

    // Create a mesh manager.
    Teuchos::RCP< MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager( 
	    mesh_blocks, getDefaultComm<int>(), mesh_dim ) );
    mesh_manager->setActiveVertices( active_vertices, 0 );

    // Partition the mesh with RCB.
    typedef RCB::zoltan_id_type zoltan_id_type;
    RCB rcb( getDefaultComm<int>(), mesh_manager, 1 );
    rcb.partition( mesh_manager->globalBoundingBox() );
    
    // Get the random numbers that were used to compute the vertex coordinates.
    std::srand( 1 );
    Teuchos::Array<double> random_numbers;
    for ( int i = 0; i < num_coords; ++i )
    {
	random_numbers.push_back( (double) std::rand() / RAND_MAX );
    }

    // Check that these are in fact the random numbers used for the vertices.
    Teuchos::ArrayRCP<const double> vertex_coords =
	mesh_blocks[0]->vertexCoordinates();
    for ( int i = 0; i < num_vertices; ++i )
    {
	TEST_ASSERT( vertex_coords[ i ] == 
		     random_numbers[mesh_dim*i] );
    }

    // Check import parameters.
    int num_import = rcb.getNumImport();

    Teuchos::ArrayView<zoltan_id_type> import_global_ids = 
	rcb.getImportGlobalIds();
    TEST_ASSERT( import_global_ids.size() == num_import );

    Teuchos::ArrayView<zoltan_id_type> import_local_ids = 
	rcb.getImportLocalIds();
    TEST_ASSERT( import_local_ids.size() == num_import );

    Teuchos::ArrayView<int> import_procs = rcb.getImportProcs();
    TEST_ASSERT( import_procs.size() == num_import );

    Teuchos::ArrayView<int> import_parts = rcb.getImportParts();
    TEST_ASSERT( import_parts.size() == num_import );

    for ( int i = 0; i < num_import; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( import_procs[i] != my_rank &&
		     import_procs[i] >= 0 &&
		     import_procs[i] < my_size );

	TEST_ASSERT( import_parts[i] == my_rank );
    }

    // Check export parameters.
    int num_export = rcb.getNumExport();

    Teuchos::ArrayView<zoltan_id_type> export_global_ids = 
	rcb.getExportGlobalIds();
    TEST_ASSERT( export_global_ids.size() == num_export );

    Teuchos::ArrayView<zoltan_id_type> export_local_ids = 
	rcb.getExportLocalIds();
    TEST_ASSERT( export_local_ids.size() == num_export );

    Teuchos::ArrayView<int> export_procs = rcb.getExportProcs();
    TEST_ASSERT( export_procs.size() == num_export );
    
    Teuchos::ArrayView<int> export_parts = rcb.getExportParts();
    TEST_ASSERT( export_parts.size() == num_export );

    for ( int i = 0; i < num_export; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( export_procs[i] == export_parts[i] );

	TEST_ASSERT( export_procs[i] != my_rank &&
		     export_procs[i] >= 0 &&
		     export_procs[i] < my_size );

	TEST_ASSERT( export_parts[i] != my_rank &&
		     export_parts[i] >= 0 &&
		     export_parts[i] < my_size );
    }

    // Check the destination proc point search.
    Teuchos::Array<double> point_0(1);
    point_0[0] = 2.0;
    TEST_ASSERT( rcb.getPointDestinationProc( point_0() ) == my_size-1 );

    Teuchos::Array<double> point_1(1);
    point_1[0] = -2.0;
    TEST_ASSERT( rcb.getPointDestinationProc( point_1() ) == 0 );

    Teuchos::Array<double> point_2(1);
    point_2[0] = 0.2;
    TEST_ASSERT( rcb.getPointDestinationProc( point_2() ) == 0 );

    Teuchos::Array<double> point_3(1);
    point_3[0] = 0.8;
    TEST_ASSERT( rcb.getPointDestinationProc( point_3() ) == my_size-1 );
}

//---------------------------------------------------------------------------//
// partial 2d mesh
TEUCHOS_UNIT_TEST( RCB, partial_2d_rcb_test )
{
    using namespace DataTransferKit;

    // create a mesh.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = build2dContainer();

    // Only some of the vertices will be partitioned.
    int mesh_dim = 2;
    int num_vertices = mesh_blocks[0]->vertexIds().size();
    int num_coords = mesh_dim * num_vertices;
    Teuchos::Array<short int> active_vertices( num_vertices, 1 );
    active_vertices[0] = 0;

    // Create a mesh manager.
    Teuchos::RCP< MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager( 
	    mesh_blocks, getDefaultComm<int>(), mesh_dim ) );
    mesh_manager->setActiveVertices( active_vertices, 0 );

    // Partition the mesh with RCB.
    typedef RCB::zoltan_id_type zoltan_id_type;
    RCB rcb( getDefaultComm<int>(), mesh_manager, 2 );
    rcb.partition( mesh_manager->globalBoundingBox() );
    
    // Get the random numbers that were used to compute the vertex coordinates.
    std::srand( 1 );
    Teuchos::Array<double> random_numbers;
    for ( int i = 0; i < num_coords; ++i )
    {
	random_numbers.push_back( (double) std::rand() / RAND_MAX );
    }

    // Check that these are in fact the random numbers used for the vertices.
    Teuchos::ArrayRCP<const double> vertex_coords =
	mesh_blocks[0]->vertexCoordinates();
    for ( int i = 0; i < num_vertices; ++i )
    {
	TEST_ASSERT( vertex_coords[ i ] == 
		     random_numbers[mesh_dim*i] );
	TEST_ASSERT( vertex_coords[ num_vertices + i ] == 
		     random_numbers[mesh_dim*i+1] );
    }

    // Get MPI parameters.
    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Check import parameters.
    int num_import = rcb.getNumImport();

    Teuchos::ArrayView<zoltan_id_type> import_global_ids = 
	rcb.getImportGlobalIds();
    TEST_ASSERT( import_global_ids.size() == num_import );

    Teuchos::ArrayView<zoltan_id_type> import_local_ids = 
	rcb.getImportLocalIds();
    TEST_ASSERT( import_local_ids.size() == num_import );

    Teuchos::ArrayView<int> import_procs = rcb.getImportProcs();
    TEST_ASSERT( import_procs.size() == num_import );

    Teuchos::ArrayView<int> import_parts = rcb.getImportParts();
    TEST_ASSERT( import_parts.size() == num_import );

    for ( int i = 0; i < num_import; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( import_procs[i] != my_rank &&
		     import_procs[i] >= 0 &&
		     import_procs[i] < my_size );

	TEST_ASSERT( import_parts[i] == my_rank );
    }

    // Check export parameters.
    int num_export = rcb.getNumExport();

    Teuchos::ArrayView<zoltan_id_type> export_global_ids = 
	rcb.getExportGlobalIds();
    TEST_ASSERT( export_global_ids.size() == num_export );

    Teuchos::ArrayView<zoltan_id_type> export_local_ids = 
	rcb.getExportLocalIds();
    TEST_ASSERT( export_local_ids.size() == num_export );

    Teuchos::ArrayView<int> export_procs = rcb.getExportProcs();
    TEST_ASSERT( export_procs.size() == num_export );
    
    Teuchos::ArrayView<int> export_parts = rcb.getExportParts();
    TEST_ASSERT( export_parts.size() == num_export );

    for ( int i = 0; i < num_export; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( export_procs[i] == export_parts[i] );

	TEST_ASSERT( export_procs[i] != my_rank &&
		     export_procs[i] >= 0 &&
		     export_procs[i] < my_size );

	TEST_ASSERT( export_parts[i] != my_rank &&
		     export_parts[i] >= 0 &&
		     export_parts[i] < my_size );
    }
}

//---------------------------------------------------------------------------//
// partial 3d mesh
TEUCHOS_UNIT_TEST( RCB, partial_3d_rcb_test )
{
    using namespace DataTransferKit;

    // create a mesh.
    Teuchos::ArrayRCP<Teuchos::RCP<MeshBlock> > mesh_blocks( 1 );
    mesh_blocks[0] = build3dContainer();

    // All of the vertices will be partitioned.
    int mesh_dim = 3;
    int num_vertices = mesh_blocks[0]->vertexIds().size();
    int num_coords = mesh_dim * num_vertices;
    Teuchos::Array<short int> active_vertices( num_vertices, 1 );
    active_vertices[0] = 0;

    // Create a mesh manager.
    Teuchos::RCP< MeshManager > mesh_manager = Teuchos::rcp( 
	new MeshManager( 
	    mesh_blocks, getDefaultComm<int>(), mesh_dim ) );
    mesh_manager->setActiveVertices( active_vertices, 0 );

    // Partition the mesh with RCB.
    typedef RCB::zoltan_id_type zoltan_id_type;
    RCB rcb( getDefaultComm<int>(), mesh_manager, 3 );
    rcb.partition( mesh_manager->globalBoundingBox() );
    
    // Get the random numbers that were used to compute the vertex coordinates.
    std::srand( 1 );
    Teuchos::Array<double> random_numbers;
    for ( int i = 0; i < num_coords; ++i )
    {
	random_numbers.push_back( (double) std::rand() / RAND_MAX );
    }

    // Check that these are in fact the random numbers used for the vertices.
    Teuchos::ArrayRCP<const double> vertex_coords =
	mesh_blocks[0]->vertexCoordinates();
    for ( int i = 0; i < num_vertices; ++i )
    {
	TEST_ASSERT( vertex_coords[ i ] == 
		     random_numbers[mesh_dim*i] );
	TEST_ASSERT( vertex_coords[ num_vertices + i ] == 
		     random_numbers[mesh_dim*i+1] );
	TEST_ASSERT( vertex_coords[ 2*num_vertices + i ] == 
		     random_numbers[mesh_dim*i+2] );
    }

    // Get MPI parameters.
    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Check import parameters.
    int num_import = rcb.getNumImport();

    Teuchos::ArrayView<zoltan_id_type> import_global_ids = 
	rcb.getImportGlobalIds();
    TEST_ASSERT( import_global_ids.size() == num_import );

    Teuchos::ArrayView<zoltan_id_type> import_local_ids = 
	rcb.getImportLocalIds();
    TEST_ASSERT( import_local_ids.size() == num_import );

    Teuchos::ArrayView<int> import_procs = rcb.getImportProcs();
    TEST_ASSERT( import_procs.size() == num_import );

    Teuchos::ArrayView<int> import_parts = rcb.getImportParts();
    TEST_ASSERT( import_parts.size() == num_import );

    for ( int i = 0; i < num_import; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( import_procs[i] != my_rank &&
		     import_procs[i] >= 0 &&
		     import_procs[i] < my_size );

	TEST_ASSERT( import_parts[i] == my_rank );
    }

    // Check export parameters.
    int num_export = rcb.getNumExport();

    Teuchos::ArrayView<zoltan_id_type> export_global_ids = 
	rcb.getExportGlobalIds();
    TEST_ASSERT( export_global_ids.size() == num_export );

    Teuchos::ArrayView<zoltan_id_type> export_local_ids = 
	rcb.getExportLocalIds();
    TEST_ASSERT( export_local_ids.size() == num_export );

    Teuchos::ArrayView<int> export_procs = rcb.getExportProcs();
    TEST_ASSERT( export_procs.size() == num_export );
    
    Teuchos::ArrayView<int> export_parts = rcb.getExportParts();
    TEST_ASSERT( export_parts.size() == num_export );

    for ( int i = 0; i < num_export; ++i )
    {
	// Check the MPI parameters.
	TEST_ASSERT( export_procs[i] == export_parts[i] );

	TEST_ASSERT( export_procs[i] != my_rank &&
		     export_procs[i] >= 0 &&
		     export_procs[i] < my_size );

	TEST_ASSERT( export_parts[i] != my_rank &&
		     export_parts[i] >= 0 &&
		     export_parts[i] < my_size );
    }
}

//---------------------------------------------------------------------------//
// end tstRCB.cpp
//---------------------------------------------------------------------------//
