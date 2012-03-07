//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tstCommIndexer.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  CommIndexer class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <Coupler_CommIndexer.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"
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
// TESTS
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( CommIndexer, duplicate_test )
{
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Communicator;

    RCP_Communicator global_comm = getDefaultComm<int>();
    RCP_Communicator local_comm = global_comm->duplicate();

    Coupler::CommIndexer<int> indexer( global_comm, local_comm );

    TEST_ASSERT( (int) indexer.size() == local_comm->getSize() );
    TEST_ASSERT( indexer.l2g( local_comm->getRank() ) == 
		 global_comm->getRank() );
}

TEUCHOS_UNIT_TEST( CommIndexer, split_test )
{
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Communicator;

    RCP_Communicator global_comm = getDefaultComm<int>();
    int rank = global_comm->getRank();
    RCP_Communicator local_comm = global_comm->split( 0, rank );

    Coupler::CommIndexer<int> indexer( global_comm, local_comm );

    TEST_ASSERT( (int) indexer.size() == local_comm->getSize() );
    TEST_ASSERT( indexer.l2g( local_comm->getRank() ) == 
		 global_comm->getRank() );
}

TEUCHOS_UNIT_TEST( CommIndexer, inverse_split_test )
{
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Communicator;

    RCP_Communicator global_comm = getDefaultComm<int>();
    int inverse_rank = global_comm->getSize() - global_comm->getRank() - 1;
    RCP_Communicator local_comm = global_comm->split( 0, inverse_rank);

    Coupler::CommIndexer<int> indexer( global_comm, local_comm );

    TEST_ASSERT( (int) indexer.size() == local_comm->getSize() );
    TEST_ASSERT( indexer.l2g( local_comm->getRank() ) == 
		 global_comm->getRank() );
}

TEUCHOS_UNIT_TEST( CommIndexer, subcommunicator_test )
{
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Communicator;

    RCP_Communicator global_comm = getDefaultComm<int>();
    std::vector<int> sub_ranks;
    for ( int n = 0; n < global_comm->getSize(); ++n )
    {
	if ( n % 2 == 0 )
	{
	    sub_ranks.push_back(n);
	}
    }
    Teuchos::ArrayView<int> sub_ranks_view( &sub_ranks[0], 
					    (int) sub_ranks.size() );
    RCP_Communicator local_comm = 
	global_comm->createSubcommunicator( sub_ranks_view );

    Coupler::CommIndexer<int> indexer( global_comm, local_comm );

    TEST_ASSERT( (int) indexer.size() == local_comm->getSize() );
    if ( global_comm->getRank() % 2 == 0 )
    {
    	TEST_ASSERT( indexer.l2g( local_comm->getRank() ) == 
		     global_comm->getRank() );
    }
    else
    {
	TEST_ASSERT( indexer.l2g( local_comm->getRank() ) == -1 );
    }
}

TEUCHOS_UNIT_TEST( CommIndexer, inverse_subcommunicator_test )
{
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Communicator;

    RCP_Communicator global_comm = getDefaultComm<int>();
    std::vector<int> inverse_sub_ranks;
    for ( int n = 0; n < global_comm->getSize(); ++n )
    {
	if ( n % 2 == 0 )
	{
	    inverse_sub_ranks.push_back(n);
	}
    }
    Teuchos::ArrayView<int> sub_ranks_view( &inverse_sub_ranks[0], 
					    (int) inverse_sub_ranks.size() );
    RCP_Communicator local_comm = 
	global_comm->createSubcommunicator( sub_ranks_view );

    Coupler::CommIndexer<int> indexer( global_comm, local_comm );

    TEST_ASSERT( (int) indexer.size() == local_comm->getSize() );
    if ( global_comm->getRank() % 2 == 0 )
    {
    	TEST_ASSERT( indexer.l2g( local_comm->getRank() ) == 
		     global_comm->getRank() );
    }
}

//---------------------------------------------------------------------------//
//                        end of tstCommIndexer.cpp
//---------------------------------------------------------------------------//
