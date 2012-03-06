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

TEUCHOS_UNIT_TEST( CommIndexer, constructor_test )
{
    Teuchos::RCP<Container> container = Teuchos::rcp( new Container() );
    Coupler::CommIndexer<Container> indexer( getDefaultComm<int>(),
					     getDefaultComm<int>(),
					     container );

    TEST_ASSERT( (int) indexer.size() == getDefaultComm<int>()->getSize() );
    TEST_ASSERT( indexer.l2g( getDefaultComm<int>()->getRank() ) ==
		 getDefaultComm<int>()->getRank() );
}


//---------------------------------------------------------------------------//
//                        end of tstCommIndexer.cpp
//---------------------------------------------------------------------------//
