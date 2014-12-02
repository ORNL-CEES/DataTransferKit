//---------------------------------------------------------------------------//
/*!
 * \file tstGeometryRCB.cpp
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

#include <DTK_Partitioner.hpp>
#include <DTK_PartitionerFactory.hpp>
#include <DTK_GeometryRCB.hpp>
#include <DTK_BoundingBox.hpp>
#include <DTK_GeometryTraits.hpp>
#include <DTK_GeometryManager.hpp>
#include <DTK_Cylinder.hpp>
#include <DTK_Box.hpp>

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
// Tests
//---------------------------------------------------------------------------//
// cylinder test
TEUCHOS_UNIT_TEST( GeometryRCB, cylinder_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Build a series of random cylinders.
    int num_cylinders = 100;
    Teuchos::ArrayRCP<Cylinder> cylinders( num_cylinders );
    Teuchos::ArrayRCP<int> gids( num_cylinders );
    for ( int i = 0; i < num_cylinders; ++i )
    {
	double length = (double) std::rand() / RAND_MAX;
	double radius = (double) std::rand() / RAND_MAX;
	double centroid_x = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	double centroid_y = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	double centroid_z = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	cylinders[i] = 
	    Cylinder( length, radius, centroid_x, centroid_y, centroid_z );
	gids[i] = i + my_rank*num_cylinders;
    }

    // Build a geometry manager.
    Teuchos::RCP<GeometryManager<Cylinder,int> > geometry_manager =
	Teuchos::rcp( new GeometryManager<Cylinder,int>( 
			  cylinders, gids, comm, 3 ) );

    Teuchos::Array<short int> active( num_cylinders, 1 );
    geometry_manager->setActiveGeometry( active );

    // Partition the geometry with GeometryRCB.
    typedef GeometryRCB<Cylinder,int>::zoltan_id_type zoltan_id_type;
    GeometryRCB<Cylinder,int> rcb( getDefaultComm<int>(), geometry_manager, 3 );
    rcb.partition( geometry_manager->globalBoundingBox() );
    
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
    Teuchos::Array<double> point_0(3);
    point_0[0] = 2.0;
    point_0[1] = 2.0;
    point_0[2] = 2.0;
    TEST_ASSERT( rcb.getPointDestinationProc( point_0() ) < my_size );

    Teuchos::Array<double> point_1(3);
    point_1[0] = -2.0;
    point_1[1] = -2.0;
    point_1[2] = -2.0;
    TEST_ASSERT( rcb.getPointDestinationProc( point_1() ) < my_size );

    Teuchos::Array<double> point_2(3);
    point_2[0] = 0.2;
    point_2[1] = 0.2;
    point_2[2] = 0.2;
    TEST_ASSERT( rcb.getPointDestinationProc( point_2() ) < my_size );

    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.8;
    point_3[1] = 0.8;
    point_3[2] = 0.8;
    TEST_ASSERT( rcb.getPointDestinationProc( point_3() ) < my_size );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( GeometryRCB, partitioner_cylinder_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Build a series of random cylinders.
    int num_cylinders = 100;
    Teuchos::ArrayRCP<Cylinder> cylinders( num_cylinders );
    Teuchos::ArrayRCP<int> gids( num_cylinders );
    for ( int i = 0; i < num_cylinders; ++i )
    {
	double length = (double) std::rand() / RAND_MAX;
	double radius = (double) std::rand() / RAND_MAX;
	double centroid_x = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	double centroid_y = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	double centroid_z = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	cylinders[i] = 
	    Cylinder( length, radius, centroid_x, centroid_y, centroid_z );
	gids[i] = i + my_rank*num_cylinders;
    }

    // Build a geometry manager.
    Teuchos::RCP<GeometryManager<Cylinder,int> > geometry_manager =
	Teuchos::rcp( new GeometryManager<Cylinder,int>( 
			  cylinders, gids, comm, 3 ) );

    Teuchos::Array<short int> active( num_cylinders, 1 );
    geometry_manager->setActiveGeometry( active );

    // Partition the geometry with GeometryRCB.
    Teuchos::RCP<Partitioner> partitioner = 
	PartitionerFactory::createGeometryPartitioner( 
	    getDefaultComm<int>(), geometry_manager, 3 );
    partitioner->partition( geometry_manager->globalBoundingBox() );
    
    // Check the destination proc point search.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 2.0;
    point_0[1] = 2.0;
    point_0[2] = 2.0;
    TEST_ASSERT( partitioner->getPointDestinationProc( point_0() ) < my_size );

    Teuchos::Array<double> point_1(3);
    point_1[0] = -2.0;
    point_1[1] = -2.0;
    point_1[2] = -2.0;
    TEST_ASSERT( partitioner->getPointDestinationProc( point_1() ) < my_size );

    Teuchos::Array<double> point_2(3);
    point_2[0] = 0.2;
    point_2[1] = 0.2;
    point_2[2] = 0.2;
    TEST_ASSERT( partitioner->getPointDestinationProc( point_2() ) < my_size );

    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.8;
    point_3[1] = 0.8;
    point_3[2] = 0.8;
    TEST_ASSERT( partitioner->getPointDestinationProc( point_3() ) < my_size );
}

//---------------------------------------------------------------------------//
// box test
TEUCHOS_UNIT_TEST( GeometryRCB, box_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Build a series of random boxes.
    int num_boxes = 100;
    Teuchos::ArrayRCP<Box> boxes( num_boxes );
    Teuchos::ArrayRCP<int> gids( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	double x_min = -(double) std::rand() / RAND_MAX + my_rank;
	double y_min = -(double) std::rand() / RAND_MAX + my_rank;
	double z_min = -(double) std::rand() / RAND_MAX + my_rank;
	double x_max =  (double) std::rand() / RAND_MAX + my_rank;
	double y_max =  (double) std::rand() / RAND_MAX + my_rank;
	double z_max =  (double) std::rand() / RAND_MAX + my_rank;
	boxes[i] = 
	    Box( x_min, y_min, z_min, x_max, y_max, z_max );
	gids[i] = i + my_rank*num_boxes;
    }

    // Build a geometry manager.
    Teuchos::RCP<GeometryManager<Box,int> > geometry_manager =
	Teuchos::rcp( new GeometryManager<Box,int>( 
			  boxes, gids, comm, 3 ) );

    Teuchos::Array<short int> active( num_boxes, 1 );
    geometry_manager->setActiveGeometry( active );

    // Partition the geometry with GeometryRCB.
    typedef GeometryRCB<Box,int>::zoltan_id_type zoltan_id_type;
    GeometryRCB<Box,int> rcb( getDefaultComm<int>(), geometry_manager, 3 );
    rcb.partition( geometry_manager->globalBoundingBox() );
    
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
    Teuchos::Array<double> point_0(3);
    point_0[0] = 2.0;
    point_0[1] = 2.0;
    point_0[2] = 2.0;
    TEST_ASSERT( rcb.getPointDestinationProc( point_0() ) < my_size );

    Teuchos::Array<double> point_1(3);
    point_1[0] = -2.0;
    point_1[1] = -2.0;
    point_1[2] = -2.0;
    TEST_ASSERT( rcb.getPointDestinationProc( point_1() ) < my_size );

    Teuchos::Array<double> point_2(3);
    point_2[0] = 0.2;
    point_2[1] = 0.2;
    point_2[2] = 0.2;
    TEST_ASSERT( rcb.getPointDestinationProc( point_2() ) < my_size );

    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.8;
    point_3[1] = 0.8;
    point_3[2] = 0.8;
    TEST_ASSERT( rcb.getPointDestinationProc( point_3() ) < my_size );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( GeometryRCB, partitioner_box_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Build a series of random boxes.
    int num_boxes = 100;
    Teuchos::ArrayRCP<Box> boxes( num_boxes );
    Teuchos::ArrayRCP<int> gids( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	double x_min = -(double) std::rand() / RAND_MAX + my_rank;
	double y_min = -(double) std::rand() / RAND_MAX + my_rank;
	double z_min = -(double) std::rand() / RAND_MAX + my_rank;
	double x_max =  (double) std::rand() / RAND_MAX + my_rank;
	double y_max =  (double) std::rand() / RAND_MAX + my_rank;
	double z_max =  (double) std::rand() / RAND_MAX + my_rank;
	boxes[i] = 
	    Box( x_min, y_min, z_min, x_max, y_max, z_max );
	gids[i] = i + my_rank*num_boxes;
    }

    // Build a geometry manager.
    Teuchos::RCP<GeometryManager<Box,int> > geometry_manager =
	Teuchos::rcp( new GeometryManager<Box,int>( 
			  boxes, gids, comm, 3 ) );

    // Partition the geometry with GeometryRCB.
    Teuchos::RCP<Partitioner> partitioner = 
	PartitionerFactory::createGeometryPartitioner( 
	    getDefaultComm<int>(), geometry_manager, 3 );
    partitioner->partition( geometry_manager->globalBoundingBox() );
    
    Teuchos::Array<short int> active( num_boxes, 1 );
    geometry_manager->setActiveGeometry( active );

    // Check the destination proc point search.
    Teuchos::Array<double> point_0(3);
    point_0[0] = 2.0;
    point_0[1] = 2.0;
    point_0[2] = 2.0;
    TEST_ASSERT( partitioner->getPointDestinationProc( point_0() ) < my_size );

    Teuchos::Array<double> point_1(3);
    point_1[0] = -2.0;
    point_1[1] = -2.0;
    point_1[2] = -2.0;
    TEST_ASSERT( partitioner->getPointDestinationProc( point_1() ) < my_size );

    Teuchos::Array<double> point_2(3);
    point_2[0] = 0.2;
    point_2[1] = 0.2;
    point_2[2] = 0.2;
    TEST_ASSERT( partitioner->getPointDestinationProc( point_2() ) < my_size );

    Teuchos::Array<double> point_3(3);
    point_3[0] = 0.8;
    point_3[1] = 0.8;
    point_3[2] = 0.8;
    TEST_ASSERT( partitioner->getPointDestinationProc( point_3() ) < my_size );
}

//---------------------------------------------------------------------------//
// callback test
TEUCHOS_UNIT_TEST( GeometryRCB, callback_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();

    // Build a series of random cylinders.
    int num_cylinders = 100;
    Teuchos::ArrayRCP<Cylinder> cylinders( num_cylinders );
    Teuchos::ArrayRCP<int> gids( num_cylinders );
    for ( int i = 0; i < num_cylinders; ++i )
    {
	double length = (double) std::rand() / RAND_MAX;
	double radius = (double) std::rand() / RAND_MAX;
	double centroid_x = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	double centroid_y = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	double centroid_z = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	cylinders[i] = 
	    Cylinder( length, radius, centroid_x, centroid_y, centroid_z );
	gids[i] = i + my_rank*num_cylinders;
    }

    // Build a geometry manager.
    Teuchos::RCP<GeometryManager<Cylinder,int> > geometry_manager =
	Teuchos::rcp( new GeometryManager<Cylinder,int>( 
			  cylinders, gids, comm, 3 ) );

    Teuchos::Array<short int> active( num_cylinders, 1 );
    active.back() = 0;
    geometry_manager->setActiveGeometry( active );
    
    Teuchos::Array<BoundingBox> boxes = geometry_manager->boundingBoxes();

    // Test the zoltan callbacks.
    void* zoltan_data = &geometry_manager;
    int error = 0;
    int num_active = num_cylinders-1;
    int num_zoltan_obj = 8*num_active;

    int num_obj = 
	GeometryRCB<Cylinder,int>::getNumberOfObjects( zoltan_data, &error );
    TEST_EQUALITY( error, 0 );
    TEST_EQUALITY( num_obj, num_zoltan_obj );

    Teuchos::Array<GeometryRCB<Cylinder,int>::zoltan_id_type> 
	global_ids( num_zoltan_obj, 0 );
    Teuchos::Array<GeometryRCB<Cylinder,int>::zoltan_id_type> 
	local_ids( num_zoltan_obj, 0 );
    int sizeGID = 1;
    int sizeLID = 1;
    GeometryRCB<Cylinder,int>::zoltan_id_ptr globalID = global_ids.getRawPtr();
    GeometryRCB<Cylinder,int>::zoltan_id_ptr localID = local_ids.getRawPtr();
    int wgt_dim = 0;
    float* obj_wgts = 0;
    GeometryRCB<Cylinder,int>::getObjectList( 
	zoltan_data, sizeGID, sizeLID, globalID, localID, wgt_dim, obj_wgts, &error );
    TEST_EQUALITY( error, 0 );
    int m = 0;
    for ( int i = 0; i < num_active; ++i )
    {
	for ( int j = 0; j < 8; ++j, ++m )
	{
	    TEST_EQUALITY( Teuchos::as<int>(local_ids[m]), m );
	    TEST_EQUALITY( Teuchos::as<int>(global_ids[m]), 8*gids[i] + j );
	}
    }

    int num_dim = 3;
    int geom_dim = GeometryRCB<Cylinder,int>::getNumGeometry( &num_dim, &error );
    TEST_EQUALITY( error, 0 );
    TEST_EQUALITY( num_dim, geom_dim );
    Teuchos::Array<double> geom_vec( num_dim*num_zoltan_obj );
    GeometryRCB<Cylinder,int>::getGeometryList( zoltan_data, sizeGID, sizeLID,
					   num_zoltan_obj, globalID, localID,
					   num_dim, geom_vec.getRawPtr(), &error );

    TEST_EQUALITY( error, 0 );
    int n = 0;
    for ( int i = 0; i < num_active; ++i )
    {
	TEST_EQUALITY( geom_vec[num_dim*n],   boxes[i].getBounds()[0] );
	TEST_EQUALITY( geom_vec[num_dim*n+1], boxes[i].getBounds()[1] );
	TEST_EQUALITY( geom_vec[num_dim*n+2], boxes[i].getBounds()[2] );
	++n;

	TEST_EQUALITY( geom_vec[num_dim*n],   boxes[i].getBounds()[3] );
	TEST_EQUALITY( geom_vec[num_dim*n+1], boxes[i].getBounds()[1] );
	TEST_EQUALITY( geom_vec[num_dim*n+2], boxes[i].getBounds()[2] );
	++n;

	TEST_EQUALITY( geom_vec[num_dim*n],   boxes[i].getBounds()[0] );
	TEST_EQUALITY( geom_vec[num_dim*n+1], boxes[i].getBounds()[4] );
	TEST_EQUALITY( geom_vec[num_dim*n+2], boxes[i].getBounds()[2] );
	++n;

	TEST_EQUALITY( geom_vec[num_dim*n],   boxes[i].getBounds()[3] );
	TEST_EQUALITY( geom_vec[num_dim*n+1], boxes[i].getBounds()[4] );
	TEST_EQUALITY( geom_vec[num_dim*n+2], boxes[i].getBounds()[2] );
	++n;

	TEST_EQUALITY( geom_vec[num_dim*n],   boxes[i].getBounds()[0] );
	TEST_EQUALITY( geom_vec[num_dim*n+1], boxes[i].getBounds()[1] );
	TEST_EQUALITY( geom_vec[num_dim*n+2], boxes[i].getBounds()[5] );
	++n;

	TEST_EQUALITY( geom_vec[num_dim*n],   boxes[i].getBounds()[3] );
	TEST_EQUALITY( geom_vec[num_dim*n+1], boxes[i].getBounds()[1] );
	TEST_EQUALITY( geom_vec[num_dim*n+2], boxes[i].getBounds()[5] );
	++n;

	TEST_EQUALITY( geom_vec[num_dim*n],   boxes[i].getBounds()[0] );
	TEST_EQUALITY( geom_vec[num_dim*n+1], boxes[i].getBounds()[4] );
	TEST_EQUALITY( geom_vec[num_dim*n+2], boxes[i].getBounds()[5] );
	++n;

	TEST_EQUALITY( geom_vec[num_dim*n],   boxes[i].getBounds()[3] );
	TEST_EQUALITY( geom_vec[num_dim*n+1], boxes[i].getBounds()[4] );
	TEST_EQUALITY( geom_vec[num_dim*n+2], boxes[i].getBounds()[5] );
	++n;
    }
}

//---------------------------------------------------------------------------//
// end tstGeometryRCB.cpp
//---------------------------------------------------------------------------//
