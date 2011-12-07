//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/test/tstPoint.cc
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  Point class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "../Point.hh"

#include "Teuchos_UnitTestHarness.hpp"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

namespace mesh {

TEUCHOS_UNIT_TEST( Point, container_test )
{
    Point<int> point(3, 4.332, 1.53, 9.87445);
    TEST_ASSERT( point.handle == 3 );
    TEST_ASSERT( point.x() == 4.332 );
    TEST_ASSERT( point.y() == 1.53 );
    TEST_ASSERT( point.z() == 9.87445 );
}

}

//---------------------------------------------------------------------------//
//                        end of tstPoint.cc
//---------------------------------------------------------------------------//
