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

#include <Coupler_Point.hpp>
#include <Coupler_SerializationTraits.hpp>

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

TEUCHOS_UNIT_TEST( Point, 1d_test )
{
    Coupler::Point<1> test_point = Coupler::point(3, 4.332);
    TEST_ASSERT( test_point.getHandle() == 3 );
    TEST_ASSERT( test_point.getCoords()[0] == 4.332 );
}

TEUCHOS_UNIT_TEST( Point, 2d_test )
{
    Coupler::Point<2> test_point = Coupler::point(3, 4.332, 1.53);
    TEST_ASSERT( test_point.getHandle() == 3 );
    TEST_ASSERT( test_point.getCoords()[0] == 4.332 );
}

TEUCHOS_UNIT_TEST( Point, 3d_test )
{
    Coupler::Point<3> test_point = Coupler::point(3, 4.332, 1.53, 9.87445);
    TEST_ASSERT( test_point.getHandle() == 3 );
    TEST_ASSERT( test_point.getCoords()[0] == 4.332 );
    TEST_ASSERT( test_point.getCoords()[1] == 1.53 );
    TEST_ASSERT( test_point.getCoords()[2] == 9.87445 );
}

TEUCHOS_UNIT_TEST( Point, 4d_test )
{
    Coupler::Point<4> test_point = Coupler::point(3, 4.332, 1.53, 9.87445, 77.3);
    TEST_ASSERT( test_point.getHandle() == 3 );
    TEST_ASSERT( test_point.getCoords()[0] == 4.332 );
    TEST_ASSERT( test_point.getCoords()[1] == 1.53 );
    TEST_ASSERT( test_point.getCoords()[2] == 9.87445 );
    TEST_ASSERT( test_point.getCoords()[3] == 77.3 );
}

TEUCHOS_UNIT_TEST( Point, 1d_serialization_test )
{
    int myRank = getDefaultComm<int>()->getRank();

    Coupler::Point<1> local_point = Coupler::point(-1, 0.0);
    
    if ( myRank == 0 )
    {
	Coupler::Point<1> broadcast_point = 
	    Coupler::point(3, 4.332);
	local_point = broadcast_point;
    }

    Teuchos::barrier<int>(*getDefaultComm<int>());

    Teuchos::broadcast<int, Coupler::Point<1> >( *getDefaultComm<int>(), 
						 0, 
						 &local_point);
    TEST_ASSERT( local_point.getHandle() == 3 );
    TEST_ASSERT( local_point.getCoords()[0] == 4.332 );
}

TEUCHOS_UNIT_TEST( Point, 2d_serialization_test )
{
    int myRank = getDefaultComm<int>()->getRank();

    Coupler::Point<2> local_point = Coupler::point(-1, 0.0, 0.0);
    
    if ( myRank == 0 )
    {
	Coupler::Point<2> broadcast_point = 
	    Coupler::point(3, 4.332, 1.53);
	local_point = broadcast_point;
    }

    Teuchos::barrier<int>(*getDefaultComm<int>());

    Teuchos::broadcast<int, Coupler::Point<2> >( *getDefaultComm<int>(), 
						 0, 
						 &local_point);
    TEST_ASSERT( local_point.getHandle() == 3 );
    TEST_ASSERT( local_point.getCoords()[0] == 4.332 );
    TEST_ASSERT( local_point.getCoords()[1] == 1.53 );
}

TEUCHOS_UNIT_TEST( Point, 3d_serialization_test )
{
    int myRank = getDefaultComm<int>()->getRank();

    Coupler::Point<3> local_point = Coupler::point(-1, 0.0, 0.0, 0.0);
    
    if ( myRank == 0 )
    {
	Coupler::Point<3> broadcast_point = 
	    Coupler::point(3, 4.332, 1.53, 9.87445);
	local_point = broadcast_point;
    }

    Teuchos::barrier<int>(*getDefaultComm<int>());

    Teuchos::broadcast<int, Coupler::Point<3> >( *getDefaultComm<int>(), 
						 0, 
						 &local_point);
    TEST_ASSERT( local_point.getHandle() == 3 );
    TEST_ASSERT( local_point.getCoords()[0] == 4.332 );
    TEST_ASSERT( local_point.getCoords()[1] == 1.53 );
    TEST_ASSERT( local_point.getCoords()[2] == 9.87445 );    
}

TEUCHOS_UNIT_TEST( Point, 4d_serialization_test )
{
    int myRank = getDefaultComm<int>()->getRank();

    Coupler::Point<4> local_point = Coupler::point(-1, 0.0, 0.0, 0.0, 0.0);
    
    if ( myRank == 0 )
    {
	Coupler::Point<4> broadcast_point = 
	    Coupler::point(3, 4.332, 1.53, 9.87445, 75.43);
	local_point = broadcast_point;
    }

    Teuchos::barrier<int>(*getDefaultComm<int>());

    Teuchos::broadcast<int, Coupler::Point<4> >( *getDefaultComm<int>(), 
						 0, 
						 &local_point);
    TEST_ASSERT( local_point.getHandle() == 3 );
    TEST_ASSERT( local_point.getCoords()[0] == 4.332 );
    TEST_ASSERT( local_point.getCoords()[1] == 1.53 );
    TEST_ASSERT( local_point.getCoords()[2] == 9.87445 );    
    TEST_ASSERT( local_point.getCoords()[3] == 75.43 );
}

//---------------------------------------------------------------------------//
//                        end of tstPoint.cpp
//---------------------------------------------------------------------------//
