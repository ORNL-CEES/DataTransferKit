//---------------------------------------------------------------------------//
/*!
 * \file tstBoundingBox.cpp
 * \author Stuart R. Slattery
 * \brief Bounding Box unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <DTK_BoundingBox.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
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
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( BoundingBox, bounding_box_test )
{
    using namespace DataTransferKit;

    // Make a bounding box.
    BoundingBox box( 3.2, -9.233, 1.3, 4.3, 0.3, 8.7 );

    // Test some points inside of it.
    double point_0[3] = { 3.7, -4, 5.4 };
    double point_1[3] = { 4.25, -7.99, 8.3 };
    double point_2[3] = { 5.4, -3, 9.4 };
    double point_3[3] = { 2.7, 0.4, 8.3 };
    TEST_ASSERT( box.pointInBox( point_0 ) );
    TEST_ASSERT( box.pointInBox( point_1 ) );
    TEST_ASSERT( !box.pointInBox( point_2 ) );
    TEST_ASSERT( !box.pointInBox( point_3 ) );
}

TEUCHOS_UNIT_TEST( BoundingBox, bounding_box_serialization_test )
{
    using namespace DataTransferKit;

    // Comm setup.
    Teuchos::RCP< const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Make a bounding box on each process.
    std::vector<BoundingBox> boxes( my_size );
    boxes[my_rank] = BoundingBox( my_rank, my_rank, my_rank,
				  my_rank+1.0, my_rank+1.0, my_rank+1.0 );

    // Reduce the boxes to an array on all processes.
    Teuchos::reduceAll<int,BoundingBox>( *comm,
					 Teuchos::REDUCE_SUM,
					 boxes.size(),
					 &boxes[0],
					 &boxes[0] );

    // Check the reduction with rank dependent coordinates.
    double rank = 0.0;
    std::vector<BoundingBox>::const_iterator box_iterator;
    for ( box_iterator = boxes.begin(); 
	  box_iterator != boxes.end(); 
	  ++box_iterator )
    {
	double point[3] = { rank+0.5, rank+0.5, rank+0.5 };
	TEST_ASSERT( box_iterator->pointInBox( point ) );

	rank += 1.0;
    }
}

//---------------------------------------------------------------------------//
// end tstBoundingBox.cpp
//---------------------------------------------------------------------------//

