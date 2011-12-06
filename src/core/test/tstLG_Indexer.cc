//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/test/tstLG_Indexer.cc
 * \author Stuart R. Slattery
 * \date   Thu Jun 16 17:00:12 2011
 * \brief  LG_Indexer unit tests
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "comm/global.hh"
#include "../LG_Indexer.hh"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// dummy application for the test
class test_app
{
  public:
    test_app() { /* ... */}
    ~test_app() {/* ... */}
};

//---------------------------------------------------------------------------//
nemesis::Communicator_t get_comm_world()
{
#ifdef COMM_MPI
    return MPI_COMM_WORLD;
#else
    return 1;
#endif
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

namespace coupler {

// Test the LG_Indexer behavior for an inverted communicator
TEUCHOS_UNIT_TEST( LG_Indexer, indexerInverted )
{
    typedef LG_Indexer        LG_Indexer_t;

    // Create the test app
    Teuchos::RCP<test_app> app_ptr( new test_app() );

    // Create the local communicator object
    nemesis::Communicator_t local_comm;

    // Split the communicator
    nemesis::split(0, nemesis::nodes()-nemesis::node()-1, local_comm);

    // Make the indexer
    LG_Indexer_t indexer(get_comm_world(), local_comm, app_ptr);

    // check the map
    TEST_ASSERT(indexer.size() == nemesis::nodes());
    for (int i = 0; i < nemesis::nodes(); ++i)
    {
	TEST_ASSERT(indexer.l2g(i) == nemesis::nodes() - i - 1);
    }
}

} // end namespace coupler

//---------------------------------------------------------------------------//
//                        end of tstLG_Indexer.cc
//---------------------------------------------------------------------------//
