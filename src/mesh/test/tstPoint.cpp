//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/test/tstPoint.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  Point class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <Mesh_Point.hpp>
#include <Mesh_SerializationTraits.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
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
// TESTS
//---------------------------------------------------------------------------//

namespace Coupler {

TEUCHOS_UNIT_TEST( Point, container_test )
{
    Point<int,double> point(3, 4.332, 1.53, 9.87445);
    TEST_ASSERT( point.handle() == 3 );
    TEST_ASSERT( point.x() == 4.332 );
    TEST_ASSERT( point.y() == 1.53 );
    TEST_ASSERT( point.z() == 9.87445 );
}

TEUCHOS_UNIT_TEST( Point, serialization_test )
{
    
    int myRank = getDefaultComm<int>()->getRank();

    Point<int,double> local_point(-1, 0.0, 0.0, 0.0);
    
    if ( myRank == 0 )
    {
	Point<int,double> broadcast_point(3, 4.332, 1.53, 9.87445);
	local_point = broadcast_point;
    }

    Teuchos::barrier<int>(*getDefaultComm<int>());

    Teuchos::broadcast<int,Point<int,double> >( *getDefaultComm<int>(), 0, &local_point);
    TEST_ASSERT( local_point.handle() == 3 );
    TEST_ASSERT( local_point.x() == 4.332 );
    TEST_ASSERT( local_point.y() == 1.53 );
    TEST_ASSERT( local_point.z() == 9.87445 );    
}

}

//---------------------------------------------------------------------------//
//                        end of tstPoint.cpp
//---------------------------------------------------------------------------//
