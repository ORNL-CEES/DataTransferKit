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

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
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
// Unit tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( GeometryManager, geometry_manager_cylinder_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_size = comm->getSize();

    // Build a series of random cylinders.
    int num_cylinders = 100;
    Teuchos::ArrayRCP<Cylinder> cylinders( num_cylinders );
    for ( int i = 0; i < num_cylinders; ++i )
    {
	double length = (double) std::rand() / RAND_MAX;
	double radius = (double) std::rand() / RAND_MAX;
	double centroid_x = (double) std::rand() / RAND_MAX - 0.5;
	double centroid_y = (double) std::rand() / RAND_MAX - 0.5;
	double centroid_z = (double) std::rand() / RAND_MAX - 0.5;
	cylinders[i] = 
	    Cylinder( length, radius, centroid_x, centroid_y, centroid_z );
    }

    // Build a geometry manager.
    GeometryManager<Cylinder> geometry_manager( cylinders, comm, 3 );

    // Check the geometry manager.
    Teuchos::ArrayRCP<Cylinder> manager_geometry = geometry_manager.geometry();
    TEST_ASSERT( manager_geometry.size() == num_cylinders );
    TEST_ASSERT( GeometryTraits<Cylinder>::dim(manager_geometry[0]) == 3 );
    TEST_ASSERT( geometry_manager.dim() == 3 );
    TEST_ASSERT( geometry_manager.comm() == comm );
    TEST_ASSERT( geometry_manager.localNumGeometry() == num_cylinders );
    TEST_ASSERT( geometry_manager.globalNumGeometry() == 
		 num_cylinders*my_size );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( GeometryManager, geometry_manager_box_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_size = comm->getSize();

    // Build a series of random boxes.
    int num_boxes = 100;
    Teuchos::ArrayRCP<Box> boxes( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
	double x_min = -(double) std::rand() / RAND_MAX;
	double y_min = -(double) std::rand() / RAND_MAX;
	double z_min = -(double) std::rand() / RAND_MAX;
	double x_max =  (double) std::rand() / RAND_MAX;
	double y_max =  (double) std::rand() / RAND_MAX;
	double z_max =  (double) std::rand() / RAND_MAX;
	boxes[i] = 
	    Box( x_min, y_min, z_min, x_max, y_max, z_max );
    }

    // Build a geometry manager.
    GeometryManager<Box> geometry_manager( boxes, comm, 3 );

    // Check the geometry manager.
    Teuchos::ArrayRCP<Box> manager_geometry = geometry_manager.geometry();
    TEST_ASSERT( manager_geometry.size() == num_boxes );
    TEST_ASSERT( GeometryTraits<Box>::dim(manager_geometry[0]) == 3 );
    TEST_ASSERT( geometry_manager.dim() == 3 );
    TEST_ASSERT( geometry_manager.comm() == comm );
    TEST_ASSERT( geometry_manager.localNumGeometry() == num_boxes );
    TEST_ASSERT( geometry_manager.globalNumGeometry() == num_boxes*my_size );
}

//---------------------------------------------------------------------------//
// end tstTransferOperator.cpp
//---------------------------------------------------------------------------//
