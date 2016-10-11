//---------------------------------------------------------------------------//
/*!
 * \file tstGeometryManager.cpp
 * \author Stuart R. Slattery
 * \brief GeometryManager unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

#include <DTK_BoundingBox.hpp>
#include <DTK_Box.hpp>
#include <DTK_Cylinder.hpp>
#include <DTK_GeometryManager.hpp>
#include <DTK_GeometryTraits.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_as.hpp>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

template <class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal>> getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp( new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Unit tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( GeometryManager, geometry_manager_cylinder_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Build a series of random cylinders.
    int num_cylinders = 100;
    Teuchos::ArrayRCP<Cylinder> cylinders( num_cylinders );
    Teuchos::ArrayRCP<unsigned long int> gids( num_cylinders );
    for ( int i = 0; i < num_cylinders; ++i )
    {
        double length = 1.0;
        double radius = 1.0;
        double centroid_x = i + my_rank;
        double centroid_y = i + my_rank;
        double centroid_z = i + my_rank;
        cylinders[i] =
            Cylinder( length, radius, centroid_x, centroid_y, centroid_z );
        gids[i] = i;
    }

    // Build a geometry manager.
    GeometryManager<Cylinder, unsigned long int> geometry_manager(
        cylinders, gids, comm, 3 );

    // Check the geometry manager.
    Teuchos::ArrayRCP<Cylinder> manager_geometry = geometry_manager.geometry();
    TEST_EQUALITY( manager_geometry.size(), num_cylinders );
    TEST_EQUALITY( GeometryTraits<Cylinder>::dim( manager_geometry[0] ), 3 );
    TEST_EQUALITY( geometry_manager.dim(), 3 );
    TEST_EQUALITY( geometry_manager.comm(), comm );
    TEST_EQUALITY( geometry_manager.localNumGeometry(), num_cylinders );
    TEST_EQUALITY( geometry_manager.globalNumGeometry(),
                   num_cylinders * my_size );

    BoundingBox local_box = geometry_manager.localBoundingBox();
    TEST_EQUALITY( local_box.getBounds()[0], my_rank - 1.0 );
    TEST_EQUALITY( local_box.getBounds()[1], my_rank - 1.0 );
    TEST_EQUALITY( local_box.getBounds()[2], my_rank - 0.5 );
    TEST_EQUALITY( local_box.getBounds()[3], my_rank + 100.0 );
    TEST_EQUALITY( local_box.getBounds()[4], my_rank + 100.0 );
    TEST_EQUALITY( local_box.getBounds()[5], my_rank + 99.5 );

    BoundingBox global_box = geometry_manager.globalBoundingBox();
    TEST_FLOATING_EQUALITY( global_box.getBounds()[0], -1.0, 1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[1], -1.0, 1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[2], -0.5, 1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[3], my_size + 99.0,
                            1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[4], my_size + 99.0,
                            1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[5], my_size + 98.5,
                            1.0e-14 );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( GeometryManager, geometry_manager_box_test )
{
    using namespace DataTransferKit;

    // Setup communication.
    Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Build a series of boxes.
    double min_val = -( my_rank + 1.0 );
    double max_val = my_rank + 1.0;
    int num_boxes = 100;
    Teuchos::ArrayRCP<Box> boxes( num_boxes );
    Teuchos::ArrayRCP<unsigned long int> gids( num_boxes );
    for ( int i = 0; i < num_boxes; ++i )
    {
        boxes[i] = Box( min_val - i, min_val - i, min_val - i, max_val + i,
                        max_val + i, max_val + i );
        gids[i] = i;
    }

    // Build a geometry manager.
    GeometryManager<Box, unsigned long int> geometry_manager( boxes, gids, comm,
                                                              3 );

    // Check the geometry manager.
    Teuchos::ArrayRCP<Box> manager_geometry = geometry_manager.geometry();
    TEST_EQUALITY( manager_geometry.size(), num_boxes );
    TEST_EQUALITY( GeometryTraits<Box>::dim( manager_geometry[0] ), 3 );
    TEST_EQUALITY( geometry_manager.dim(), 3 );
    TEST_EQUALITY( geometry_manager.comm(), comm );
    TEST_EQUALITY( geometry_manager.localNumGeometry(), num_boxes );
    TEST_EQUALITY( geometry_manager.globalNumGeometry(), num_boxes * my_size );

    double dbl_size = my_size;

    BoundingBox local_box = geometry_manager.localBoundingBox();
    TEST_EQUALITY( local_box.getBounds()[0], min_val - 99.0 );
    TEST_EQUALITY( local_box.getBounds()[1], min_val - 99.0 );
    TEST_EQUALITY( local_box.getBounds()[2], min_val - 99.0 );
    TEST_EQUALITY( local_box.getBounds()[3], max_val + 99.0 );
    TEST_EQUALITY( local_box.getBounds()[4], max_val + 99.0 );
    TEST_EQUALITY( local_box.getBounds()[5], max_val + 99.0 );

    BoundingBox global_box = geometry_manager.globalBoundingBox();
    TEST_FLOATING_EQUALITY( global_box.getBounds()[0], -dbl_size - 99.0,
                            1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[1], -dbl_size - 99.0,
                            1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[2], -dbl_size - 99.0,
                            1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[3], dbl_size + 99.0,
                            1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[4], dbl_size + 99.0,
                            1.0e-14 );
    TEST_FLOATING_EQUALITY( global_box.getBounds()[5], dbl_size + 99.0,
                            1.0e-14 );
}

//---------------------------------------------------------------------------//
// end tstTransferOperator.cpp
//---------------------------------------------------------------------------//
