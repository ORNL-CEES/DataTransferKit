//---------------------------------------------------------------------------//
/*!
 * \file   tstCommTools.cpp
 * \author Stuart Slattery
 * \brief  CommTools class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <DTK_CommTools.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// Get the default communicator.
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
// Tests.
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( CommTools, comm_world_test )
{
    using namespace DataTransferKit;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm comm_world;
    CommTools::getMpiCommWorld( comm_world );

    RCP_Comm default_comm = getDefaultComm<int>();
    TEST_ASSERT( comm_world->getRank() == default_comm->getRank() );
    TEST_ASSERT( comm_world->getSize() == default_comm->getSize() );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommTools, equal_test_1 )
{
    using namespace DataTransferKit;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm comm_A = getDefaultComm<int>();
    RCP_Comm comm_B = comm_A->duplicate();

    TEST_ASSERT( CommTools::equal( comm_A, comm_B ) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommTools, equal_test_2 )
{
    using namespace DataTransferKit;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm comm_A = getDefaultComm<int>();
    Teuchos::Array<int> sub_ranks;
    for ( int n = 0; n < comm_A->getSize(); ++n )
    {
	if ( n % 2 == 0 )
	{
	    sub_ranks.push_back(n);
	}
    }
    RCP_Comm comm_B = 
	comm_A->createSubcommunicator( sub_ranks() );

    TEST_ASSERT( !CommTools::equal( comm_A, comm_B ) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommTools, union_test_1 )
{
    using namespace DataTransferKit;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm comm_A = getDefaultComm<int>();
    RCP_Comm comm_B = comm_A->duplicate();
    RCP_Comm comm_union;
    CommTools::unite( comm_A, comm_B, comm_union );

    TEST_ASSERT( CommTools::equal( comm_A, comm_union ) );
    TEST_ASSERT( CommTools::equal( comm_B, comm_union ) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommTools, union_test_2 )
{
    using namespace DataTransferKit;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm comm_A = getDefaultComm<int>();
    Teuchos::Array<int> sub_ranks;
    for ( int n = 0; n < comm_A->getSize(); ++n )
    {
	if ( n % 2 == 0 )
	{
	    sub_ranks.push_back(n);
	}
    }
    RCP_Comm comm_B = 
	comm_A->createSubcommunicator( sub_ranks() );

    RCP_Comm comm_union;
    CommTools::unite( comm_A, comm_B, comm_union );
    TEST_ASSERT( CommTools::equal( comm_A, comm_union ) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommTools, intersect_test_1 )
{
    using namespace DataTransferKit;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm comm_A = getDefaultComm<int>();
    RCP_Comm comm_B = comm_A->duplicate();
    RCP_Comm comm_intersect;
    CommTools::intersect( comm_A, comm_B, comm_intersect );

    TEST_ASSERT( CommTools::equal( comm_A, comm_intersect ) );
    TEST_ASSERT( CommTools::equal( comm_B, comm_intersect ) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommTools, intersect_test_2 )
{
    using namespace DataTransferKit;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm comm_A = getDefaultComm<int>();
    Teuchos::Array<int> sub_ranks;
    for ( int n = 0; n < comm_A->getSize(); ++n )
    {
	if ( n % 2 == 0 )
	{
	    sub_ranks.push_back(n);
	}
    }
    RCP_Comm comm_B = 
	comm_A->createSubcommunicator( sub_ranks() );

    RCP_Comm comm_intersect;
    CommTools::intersect( comm_A, comm_B, comm_intersect );
    TEST_ASSERT( CommTools::equal( comm_B, comm_intersect ) );
}

//---------------------------------------------------------------------------//
// end tstCommTools.cpp
//---------------------------------------------------------------------------//
