//---------------------------------------------------------------------------//
/*!
 * \file tstGeometryLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief Bounding Box unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_Box.hpp>
#include <DTK_BasicGeometryLocalMap.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>

//---------------------------------------------------------------------------//
// Global test variables.
//---------------------------------------------------------------------------//
int num_rand = 1000;

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( BasicGeometryLocalMap, local_map_test )
{
    using namespace DataTransferKit;

    // make a box.
    double x_min = 3.2;
    double y_min = -9.233;
    double z_min = 1.3;
    double x_max = 4.3;
    double y_max = 0.3;
    double z_max = 8.7;
    Entity box = Box(  0, 0, 0, x_min, y_min, z_min, x_max, y_max, z_max );

    // Make a local map.
    Teuchos::RCP<EntityLocalMap> local_map = 
	Teuchos::rcp( new BasicGeometryLocalMap() );

    // Compute the measure.
    TEST_FLOATING_EQUALITY( local_map->measure(box), 77.5986, 1.0e-4 );

    // Test some random points inside of it.
    Teuchos::Tuple<double,6> box_bounds;
    box.boundingBox( box_bounds );
    Teuchos::Array<double> point(3);
    Teuchos::Array<double> ref_point(3);
    Teuchos::ParameterList plist;
    bool point_inclusion = false;
    bool safe_to_map = false;
    bool map_ok = false;
    for ( int i = 0; i < num_rand; ++i )
    {
	point[0] = 2.0 * (double) std::rand() / RAND_MAX + 3.0;
	point[1] = 12.0 * (double) std::rand() / RAND_MAX - 11.0;
	point[2] = 9.0 * (double) std::rand() / RAND_MAX;

	safe_to_map = local_map->isSafeToMapToReferenceFrame( box, point() );
	map_ok = local_map->mapToReferenceFrame( box, point, ref_point() );
	TEST_ASSERT( map_ok );
	point_inclusion = local_map->checkPointInclusion(box, ref_point() );
	if ( box_bounds[0] <= ref_point[0] &&
	     box_bounds[1] <= ref_point[1] &&
	     box_bounds[2] <= ref_point[2] &&
	     box_bounds[3] >= ref_point[0] &&
	     box_bounds[4] >= ref_point[1] &&
	     box_bounds[5] >= ref_point[2] )
	{
	    TEST_ASSERT( point_inclusion );
	    TEST_ASSERT( safe_to_map );
	}
	else
	{
	    TEST_ASSERT( !point_inclusion );
	    TEST_ASSERT( !safe_to_map );
	}
    }
}

//---------------------------------------------------------------------------//
// end tstGeometryLocalMap.cpp
//---------------------------------------------------------------------------//

