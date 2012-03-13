//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/test/tstBoundingBox.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  BoundingBox class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <DataTransferKit_Point.hpp>
#include <DataTransferKit_BoundingBox.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Tuple.hpp"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

namespace DataTransferKit {

TEUCHOS_UNIT_TEST( Point, container_test )
{
    BoundingBox<int,double> box(1.3, 5.44, -0.32, 98.4, 5.4, 10.2);
    TEST_ASSERT( box.domain()[0] == 1.3 );
    TEST_ASSERT( box.domain()[1] == 5.44 );
    TEST_ASSERT( box.domain()[2] == -0.32 );
    TEST_ASSERT( box.domain()[3] == 98.4 );
    TEST_ASSERT( box.domain()[4] == 5.4 );
    TEST_ASSERT( box.domain()[5] == 10.2 );

    Point<int,double> inside_point(3, 4.332, 1.53, 9.87445);
    TEST_ASSERT( box.point_query(inside_point) );

    Point<int,double> outside_point(6, -43.3, 153.3, 3.2);
    TEST_ASSERT( !box.point_query(outside_point) );
}

}

//---------------------------------------------------------------------------//
//                        end of tstBoundingBox.cpp
//---------------------------------------------------------------------------//
