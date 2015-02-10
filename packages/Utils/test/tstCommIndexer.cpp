//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
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

#include <DTK_CommIndexer.hpp>

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
    using namespace DataTransferKit;

    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm global_comm = getDefaultComm<int>();
    RCP_Comm local_comm = global_comm->duplicate();

    CommIndexer indexer( global_comm, local_comm );

    TEST_ASSERT( (int) indexer.size() == local_comm->getSize() );
    TEST_ASSERT( indexer.l2g( local_comm->getRank() ) == 
		 global_comm->getRank() );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommIndexer, split_test )
{
    using namespace DataTransferKit;

    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm global_comm = getDefaultComm<int>();
    int rank = global_comm->getRank();
    RCP_Comm local_comm = global_comm->split( 0, rank );

    CommIndexer indexer( global_comm, local_comm );

    TEST_ASSERT( (int) indexer.size() == local_comm->getSize() );
    TEST_ASSERT( indexer.l2g( local_comm->getRank() ) == 
		 global_comm->getRank() );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommIndexer, inverse_split_test )
{
    using namespace DataTransferKit;

    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm global_comm = getDefaultComm<int>();
    int inverse_rank = global_comm->getSize() - global_comm->getRank() - 1;
    RCP_Comm local_comm = global_comm->split( 0, inverse_rank);

    CommIndexer indexer( global_comm, local_comm );

    TEST_ASSERT( (int) indexer.size() == local_comm->getSize() );
    TEST_ASSERT( indexer.l2g( local_comm->getRank() ) == 
		 global_comm->getRank() );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommIndexer, subcommunicator_test )
{
    using namespace DataTransferKit;

    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm global_comm = getDefaultComm<int>();
    Teuchos::Array<int> sub_ranks;
    for ( int n = 0; n < global_comm->getSize(); ++n )
    {
	if ( n % 2 == 0 )
	{
	    sub_ranks.push_back(n);
	}
    }
    RCP_Comm local_comm = 
	global_comm->createSubcommunicator( sub_ranks() );

    CommIndexer indexer( global_comm, local_comm );

    if ( global_comm->getRank() % 2 == 0 )
    {
	TEST_ASSERT( indexer.isValid() );
	TEST_ASSERT( Teuchos::nonnull(local_comm) );
    	TEST_EQUALITY( (int) indexer.size(), local_comm->getSize() );
    	TEST_EQUALITY( indexer.l2g( local_comm->getRank() ), 
    		       global_comm->getRank() );
    }
    else
    {
	TEST_ASSERT( !indexer.isValid() );
	TEST_ASSERT( Teuchos::is_null(local_comm) );
    }
    TEST_EQUALITY( indexer.l2g( -32 ), -1 );
}

//---------------------------------------------------------------------------//
//                        end of tstCommIndexer.cpp
//---------------------------------------------------------------------------//
