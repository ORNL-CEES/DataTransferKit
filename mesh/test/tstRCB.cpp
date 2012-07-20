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
#include <DTK_MeshTraits.hpp>
#include <DTK_MeshManager.hpp>
#include <DTK_MeshContainer.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
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
// Mesh container creation functions.
//---------------------------------------------------------------------------//
// 1d mesh
DataTransferKit::MeshContainer<int> build1dContainer()
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Make some random nodes.
    int node_dim = 1;
    int num_rand = node_dim*my_size;
    std::srand( 1 );
    Teuchos::ArrayRCP<int> node_handles( my_size );
    Teuchos::ArrayRCP<double> coords( num_rand );
    for ( int i = 0; i < my_size; ++i )
    {
	node_handles[i] = i*my_size + my_rank;
	coords[ i ] = (double) std::rand() / RAND_MAX;
    }

    // Empty element vectors. We only need nodes for these tests.
    Teuchos::ArrayRCP<int> element_handles;
    Teuchos::ArrayRCP<int> element_connectivity;
    Teuchos::ArrayRCP<std::size_t> permutation_list;

    return MeshContainer<int>( node_dim, node_handles, coords, 
			       DTK_LINE_SEGMENT, 2,
			       element_handles, element_connectivity,
			       permutation_list );
}

//---------------------------------------------------------------------------//
// 2d mesh
DataTransferKit::MeshContainer<int> build2dContainer()
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Make some random nodes.
    int node_dim = 2;
    int num_rand = node_dim*my_size;
    std::srand( 1 );
    Teuchos::ArrayRCP<int> node_handles( my_size );
    Teuchos::ArrayRCP<double> coords( num_rand );
    for ( int i = 0; i < my_size; ++i )
    {
	node_handles[i] = i*my_size + my_rank;
	coords[ i ] = (double) std::rand() / RAND_MAX;
	coords[ my_size + i ] = (double) std::rand() / RAND_MAX;
    }

    // Empty element vectors. We only need nodes for these tests.
    Teuchos::ArrayRCP<int> element_handles;
    Teuchos::ArrayRCP<int> element_connectivity;
    Teuchos::ArrayRCP<std::size_t> permutation_list;

    return MeshContainer<int>( node_dim, node_handles, coords, 
			       DTK_TRIANGLE, 3,
			       element_handles, element_connectivity,
			       permutation_list );
}

//---------------------------------------------------------------------------//
// 3d mesh
DataTransferKit::MeshContainer<int> build3dContainer()
{
    using namespace DataTransferKit;

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    // Make some random nodes.
    int node_dim = 3;
    int num_rand = node_dim*my_size;
    std::srand( 1 );
    Teuchos::ArrayRCP<int> node_handles( my_size );
    Teuchos::ArrayRCP<double> coords( num_rand );
    for ( int i = 0; i < my_size; ++i )
    {
	node_handles[i] = i*my_size + my_rank;
	coords[ i ] = (double) std::rand() / RAND_MAX;
	coords[ my_size + i ] = (double) std::rand() / RAND_MAX;
	coords[ 2*my_size + i ] = (double) std::rand() / RAND_MAX;
    }

    // Empty element vectors. We only need nodes for these tests.
    Teuchos::ArrayRCP<int> element_handles;
    Teuchos::ArrayRCP<int> element_connectivity;
    Teuchos::ArrayRCP<std::size_t> permutation_list;

    return MeshContainer<int>( node_dim, node_handles, coords, 
			       DTK_HEXAHEDRON, 8,
			       element_handles, element_connectivity,
			       permutation_list );
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
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP< MeshType > mesh_blocks( 1 );
    mesh_blocks[0] = build1dContainer();

    // All of the nodes will be partitioned.
    int mesh_dim = 1;
    int num_nodes = Tools::numNodes( mesh_blocks[0] );
    int num_coords = mesh_dim * num_nodes;
    Teuchos::Array<short int> active_nodes( num_nodes, 1 );

    // Create a mesh manager.
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>( 
	    mesh_blocks, getDefaultComm<int>(), mesh_dim ) );
    mesh_manager->setActiveNodes( active_nodes, 0 );

    // Partition the mesh with RCB.
    typedef RCB<MeshType>::zoltan_id_type zoltan_id_type;
    RCB<MeshType> rcb( mesh_manager );
    rcb.partition();
    
    // Get the random numbers that were used to compute the node coordinates.
    std::srand( 1 );
    Teuchos::Array<double> random_numbers;
    for ( int i = 0; i < num_coords; ++i )
    {
	random_numbers.push_back( (double) std::rand() / RAND_MAX );
    }

    // Check that these are in fact the random numbers used for the nodes.
    typename MT::const_coordinate_iterator coord_iterator 
	= MT::coordsBegin( mesh_blocks[0] );
    for ( int i = 0; i < num_nodes; ++i )
    {
	TEST_ASSERT( coord_iterator[ i ] == 
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
}

//---------------------------------------------------------------------------//
// 2d mesh
TEUCHOS_UNIT_TEST( RCB, 2d_rcb_test )
{
    using namespace DataTransferKit;

    // create a mesh.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP< MeshType > mesh_blocks( 1 );
    mesh_blocks[0] = build2dContainer();

    // All of the nodes will be partitioned.
    int mesh_dim = 2;
    int num_nodes = Tools::numNodes( mesh_blocks[0] );
    int num_coords = mesh_dim * num_nodes;
    Teuchos::Array<short int> active_nodes( num_nodes, 1 );

    // Create a mesh manager.
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>( 
	    mesh_blocks, getDefaultComm<int>(), mesh_dim ) );
    mesh_manager->setActiveNodes( active_nodes, 0 );

    // Partition the mesh with RCB.
    typedef RCB<MeshType>::zoltan_id_type zoltan_id_type;
    RCB<MeshType> rcb( mesh_manager );
    rcb.partition();
    
    // Get the random numbers that were used to compute the node coordinates.
    std::srand( 1 );
    Teuchos::Array<double> random_numbers;
    for ( int i = 0; i < num_coords; ++i )
    {
	random_numbers.push_back( (double) std::rand() / RAND_MAX );
    }

    // Check that these are in fact the random numbers used for the nodes.
    typename MT::const_coordinate_iterator coord_iterator 
	= MT::coordsBegin( mesh_blocks[0] );
    for ( int i = 0; i < num_nodes; ++i )
    {
	TEST_ASSERT( coord_iterator[ i ] == 
		     random_numbers[mesh_dim*i] );
	TEST_ASSERT( coord_iterator[ num_nodes + i ] == 
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
// 3d mesh
TEUCHOS_UNIT_TEST( RCB, 3d_rcb_test )
{
    using namespace DataTransferKit;

    // create a mesh.
    typedef MeshContainer<int> MeshType;
    typedef MeshTraits< MeshType > MT;
    typedef MeshTools< MeshType > Tools;
    Teuchos::ArrayRCP< MeshType > mesh_blocks( 1 );
    mesh_blocks[0] = build3dContainer();

    // All of the nodes will be partitioned.
    int mesh_dim = 3;
    int num_nodes = Tools::numNodes( mesh_blocks[0] );
    int num_coords = mesh_dim * num_nodes;
    Teuchos::Array<short int> active_nodes( num_nodes, 1 );

    // Create a mesh manager.
    Teuchos::RCP< MeshManager<MeshType> > mesh_manager = Teuchos::rcp( 
	new MeshManager<MeshType>( 
	    mesh_blocks, getDefaultComm<int>(), mesh_dim ) );
    mesh_manager->setActiveNodes( active_nodes, 0 );

    // Partition the mesh with RCB.
    typedef RCB<MeshType>::zoltan_id_type zoltan_id_type;
    RCB<MeshType> rcb( mesh_manager );
    rcb.partition();
    
    // Get the random numbers that were used to compute the node coordinates.
    std::srand( 1 );
    Teuchos::Array<double> random_numbers;
    for ( int i = 0; i < num_coords; ++i )
    {
	random_numbers.push_back( (double) std::rand() / RAND_MAX );
    }

    // Check that these are in fact the random numbers used for the nodes.
    typename MT::const_coordinate_iterator coord_iterator 
	= MT::coordsBegin( mesh_blocks[0] );
    for ( int i = 0; i < num_nodes; ++i )
    {
	TEST_ASSERT( coord_iterator[ i ] == 
		     random_numbers[mesh_dim*i] );
	TEST_ASSERT( coord_iterator[ num_nodes + i ] == 
		     random_numbers[mesh_dim*i+1] );
	TEST_ASSERT( coord_iterator[ 2*num_nodes + i ] == 
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
