//---------------------------------------------------------------------------//
/*!
 * \file tstGeometryManager.cpp
 * \author Stuart R. Slattery
 * \brief GeometryManager unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_GeometryManager.hpp>
#include <DTK_GeometryTraits.hpp>
#include <DTK_Cylinder.hpp>
#include <DTK_Box.hpp>
#include <DTK_BoundingBox.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_as.hpp>

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
// Unit tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( GeometryManager, geometry_manager_cylinder_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Build a series of random cylinders.
    int num_cylinders = 100;
    Teuchos::ArrayRCP<Cylinder> cylinders( num_cylinders );
    Teuchos::ArrayRCP<unsigned long int> gids( num_cylinders );
    for ( int i = 0; i < num_cylinders; ++i )
    {
	double length = (double) std::rand() / RAND_MAX;
	double radius = (double) std::rand() / RAND_MAX;
	double centroid_x = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	double centroid_y = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	double centroid_z = (double) std::rand() / RAND_MAX - 0.5 + my_rank;
	cylinders[i] = 
	    Cylinder( length, radius, centroid_x, centroid_y, centroid_z );
	gids[i] = i;
    }

    // Build a geometry manager.
    GeometryManager<Cylinder,unsigned long int> geometry_manager( cylinders, gids, comm, 3 );

    // Check the geometry manager.
    Teuchos::ArrayRCP<Cylinder> manager_geometry = geometry_manager.geometry();
    TEST_ASSERT( manager_geometry.size() == num_cylinders );
    TEST_ASSERT( GeometryTraits<Cylinder>::dim(manager_geometry[0]) == 3 );
    TEST_ASSERT( geometry_manager.dim() == 3 );
    TEST_ASSERT( geometry_manager.comm() == comm );
    TEST_ASSERT( geometry_manager.localNumGeometry() == num_cylinders );
    TEST_ASSERT( geometry_manager.globalNumGeometry() == 
		 num_cylinders*my_size );

    // Check the bounding box functions.
    double x_min = Teuchos::ScalarTraits<double>::rmax();
    double y_min = Teuchos::ScalarTraits<double>::rmax();
    double z_min = Teuchos::ScalarTraits<double>::rmax();
    double x_max = -Teuchos::ScalarTraits<double>::rmax();
    double y_max = -Teuchos::ScalarTraits<double>::rmax();
    double z_max = -Teuchos::ScalarTraits<double>::rmax();

    Teuchos::Array<BoundingBox> cylinder_boxes = 
	geometry_manager.boundingBoxes();
    Teuchos::ArrayRCP<unsigned long int> cylinder_gids = geometry_manager.gids();
    for ( int i = 0; i < num_cylinders; ++i )
    {
	TEST_ASSERT( cylinder_gids[i] == gids[i] );

	for ( int d = 0; d < 6; ++d )
	{
	    TEST_ASSERT( cylinder_boxes[i].getBounds()[d] ==
			 cylinders[i].boundingBox().getBounds()[d] );
	}

	if ( cylinders[i].boundingBox().getBounds()[0] < x_min )
	    x_min = cylinders[i].boundingBox().getBounds()[0];

	if ( cylinders[i].boundingBox().getBounds()[1] < y_min )
	    y_min = cylinders[i].boundingBox().getBounds()[1];

	if ( cylinders[i].boundingBox().getBounds()[2] < z_min )
	    z_min = cylinders[i].boundingBox().getBounds()[2];

	if ( cylinders[i].boundingBox().getBounds()[3] > x_max )
	    x_max = cylinders[i].boundingBox().getBounds()[3];

	if ( cylinders[i].boundingBox().getBounds()[4] > y_max )
	    y_max = cylinders[i].boundingBox().getBounds()[4];

	if ( cylinders[i].boundingBox().getBounds()[5] > z_max )
	    z_max = cylinders[i].boundingBox().getBounds()[5];
    }

    BoundingBox local_box = geometry_manager.localBoundingBox();
    TEST_ASSERT( local_box.getBounds()[0] == x_min );
    TEST_ASSERT( local_box.getBounds()[1] == y_min );
    TEST_ASSERT( local_box.getBounds()[2] == z_min );
    TEST_ASSERT( local_box.getBounds()[3] == x_max );
    TEST_ASSERT( local_box.getBounds()[4] == y_max );
    TEST_ASSERT( local_box.getBounds()[5] == z_max );

    BoundingBox global_box = geometry_manager.globalBoundingBox();
    TEST_FLOATING_EQUALITY( global_box.getBounds()[0] , x_min-my_rank, 1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[1] , y_min-my_rank, 1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[2] , z_min-my_rank, 1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[3] , x_max+my_size-my_rank-1, 1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[4] , y_max+my_size-my_rank-1, 1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[5] , z_max+my_size-my_rank-1, 1.0e-14 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( GeometryManager, geometry_manager_box_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Build a series of random boxes.
    int num_boxes = 100;
    Teuchos::ArrayRCP<Box> boxes( num_boxes );
    Teuchos::ArrayRCP<unsigned long int> gids( num_boxes );
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
	gids[i] = i;
    }

    // Build a geometry manager.
    GeometryManager<Box,unsigned long int> geometry_manager( boxes, gids, comm, 3 );

    // Check the geometry manager.
    Teuchos::ArrayRCP<Box> manager_geometry = geometry_manager.geometry();
    TEST_ASSERT( manager_geometry.size() == num_boxes );
    TEST_ASSERT( GeometryTraits<Box>::dim(manager_geometry[0]) == 3 );
    TEST_ASSERT( geometry_manager.dim() == 3 );
    TEST_ASSERT( geometry_manager.comm() == comm );
    TEST_ASSERT( geometry_manager.localNumGeometry() == num_boxes );
    TEST_ASSERT( geometry_manager.globalNumGeometry() == num_boxes*my_size );

    // Check the bounding box functions.
    double x_min = Teuchos::ScalarTraits<double>::rmax();
    double y_min = Teuchos::ScalarTraits<double>::rmax();
    double z_min = Teuchos::ScalarTraits<double>::rmax();
    double x_max = -Teuchos::ScalarTraits<double>::rmax();
    double y_max = -Teuchos::ScalarTraits<double>::rmax();
    double z_max = -Teuchos::ScalarTraits<double>::rmax();

    Teuchos::Array<BoundingBox> box_boxes = 
	geometry_manager.boundingBoxes();
    Teuchos::ArrayRCP<unsigned long int> box_gids = geometry_manager.gids();
    for ( int i = 0; i < num_boxes; ++i )
    {
	TEST_ASSERT( box_gids[i] == gids[i] );

	for ( int d = 0; d < 6; ++d )
	{
	    TEST_ASSERT( box_boxes[i].getBounds()[d] ==
			 boxes[i].boundingBox().getBounds()[d] );
	}

	if ( boxes[i].boundingBox().getBounds()[0] < x_min )
	    x_min = boxes[i].boundingBox().getBounds()[0];

	if ( boxes[i].boundingBox().getBounds()[1] < y_min )
	    y_min = boxes[i].boundingBox().getBounds()[1];

	if ( boxes[i].boundingBox().getBounds()[2] < z_min )
	    z_min = boxes[i].boundingBox().getBounds()[2];

	if ( boxes[i].boundingBox().getBounds()[3] > x_max )
	    x_max = boxes[i].boundingBox().getBounds()[3];

	if ( boxes[i].boundingBox().getBounds()[4] > y_max )
	    y_max = boxes[i].boundingBox().getBounds()[4];

	if ( boxes[i].boundingBox().getBounds()[5] > z_max )
	    z_max = boxes[i].boundingBox().getBounds()[5];
    }

    BoundingBox local_box = geometry_manager.localBoundingBox();
    TEST_ASSERT( local_box.getBounds()[0] == x_min );
    TEST_ASSERT( local_box.getBounds()[1] == y_min );
    TEST_ASSERT( local_box.getBounds()[2] == z_min );
    TEST_ASSERT( local_box.getBounds()[3] == x_max );
    TEST_ASSERT( local_box.getBounds()[4] == y_max );
    TEST_ASSERT( local_box.getBounds()[5] == z_max );

    BoundingBox global_box = geometry_manager.globalBoundingBox();
    TEST_FLOATING_EQUALITY( global_box.getBounds()[0] , x_min-my_rank, 1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[1] , y_min-my_rank, 1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[2] , z_min-my_rank, 1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[3] , x_max+my_size-my_rank-1, 1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[4] , y_max+my_size-my_rank-1, 1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[5] , z_max+my_size-my_rank-1, 1.0e-14 );
}

//---------------------------------------------------------------------------//
// end tstTransferOperator.cpp
//---------------------------------------------------------------------------//
