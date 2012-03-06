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

// A simple container to test the indexer with.
class Container
{
  public:

    Container()
    { /* ... */ }

    ~Container()
    { /* ... */ }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( CommIndexer, duplicate_test )
{
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Communicator;

    RCP_Communicator comm_1 = getDefaultComm<int>();
    RCP_Communicator comm_2 = comm_1->duplicate();

    Teuchos::RCP<Container> container = Teuchos::rcp( new Container() );
    Coupler::CommIndexer<Container> indexer( comm_1, comm_2, container );

    TEST_ASSERT( (int) indexer.size() == comm_2->getSize() );
    TEST_ASSERT( indexer.l2g( comm_2->getRank() ) == comm_1->getRank() );
}

TEUCHOS_UNIT_TEST( CommIndexer, inverse_duplicate_test )
{
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Communicator;

    RCP_Communicator comm_1 = getDefaultComm<int>();
    int inverse_rank = comm_1->getSize() - comm_1->getRank() - 1;
    RCP_Communicator comm_2 = comm_1->split( 0, inverse_rank);

    Teuchos::RCP<Container> container = Teuchos::rcp( new Container() );
    Coupler::CommIndexer<Container> indexer( comm_1, comm_2, container );

    TEST_ASSERT( (int) indexer.size() == comm_2->getSize() );
    TEST_ASSERT( indexer.l2g( comm_2->getRank() ) == comm_1->getRank() );
}

TEUCHOS_UNIT_TEST( CommIndexer, subcommunicator_test )
{
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Communicator;

    RCP_Communicator comm_1 = getDefaultComm<int>();
    std::vector<int> sub_ranks;
    for ( int n = 0; n < comm_1->getSize(); ++n )
    {
	if ( n % 2 == 0 )
	{
	    sub_ranks.push_back(n);
	}
    }
    Teuchos::ArrayView<int> sub_ranks_view( &sub_ranks[0], 
					    (int) sub_ranks.size() );
    RCP_Communicator comm_2 = comm_1->createSubcommunicator( sub_ranks_view );

    Teuchos::RCP<Container> container = Teuchos::rcp( new Container() );
    Coupler::CommIndexer<Container> indexer( comm_1, comm_2, container );

    TEST_ASSERT( (int) indexer.size() == comm_2->getSize() );
    if ( comm_1->getRank() % 2 == 0 )
    {
	TEST_ASSERT( indexer.l2g( comm_2->getRank() ) == comm_1->getRank() );
    }
}

//---------------------------------------------------------------------------//
//                        end of tstCommIndexer.cpp
//---------------------------------------------------------------------------//
