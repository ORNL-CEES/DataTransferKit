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
    CommTools::getCommWorld( comm_world );

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

    if ( comm_A->getSize() > 1 )
    {
        TEST_ASSERT( !CommTools::equal( comm_A, comm_B ) );
    }
    else
    {
        TEST_ASSERT( CommTools::equal( comm_A, comm_B ) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommTools, equal_test_3 )
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

    RCP_Comm comm_global = comm_A;

    if ( comm_A->getSize() > 1 )
    {
        TEST_ASSERT( !CommTools::equal( comm_A, comm_B, comm_global ) );
    }
    else
    {
        TEST_ASSERT( CommTools::equal( comm_A, comm_B, comm_global ) );
    }
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

    if ( comm_A->getSize() > 1 )
    {
        TEST_ASSERT( !CommTools::equal( comm_union, comm_B ) );
    }
    else
    {
        TEST_ASSERT( CommTools::equal( comm_union, comm_B ) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommTools, union_test_3 )
{
    using namespace DataTransferKit;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm comm_A = getDefaultComm<int>();
    int inverse_rank = comm_A->getSize() - comm_A->getRank() - 1;
    RCP_Comm comm_B = comm_A->split( 0, inverse_rank);

    RCP_Comm comm_union;
    CommTools::unite( comm_A, comm_B, comm_union );
    TEST_ASSERT( CommTools::equal( comm_A, comm_union ) );
    TEST_ASSERT( CommTools::equal( comm_B, comm_union ) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommTools, union_test_4 )
{
    using namespace DataTransferKit;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm comm_default = getDefaultComm<int>();
    Teuchos::Array<int> sub_ranks_a, sub_ranks_b;
    for ( int n = 0; n < comm_default->getSize(); ++n )
    {
        if ( n % 2 == 0 )
        {
            sub_ranks_a.push_back(n);
        }
        else
        {
            sub_ranks_b.push_back(n);
        }
    }

    RCP_Comm comm_A =
        comm_default->createSubcommunicator( sub_ranks_a() );

    RCP_Comm comm_B =
        comm_default->createSubcommunicator( sub_ranks_b() );

    RCP_Comm comm_union;
    CommTools::unite( comm_A, comm_B, comm_union );

    TEST_ASSERT( CommTools::equal( comm_default, comm_union ) );

    if ( comm_default->getSize() > 1 )
    {
        TEST_ASSERT( !CommTools::equal( comm_A, comm_union ) );
        TEST_ASSERT( !CommTools::equal( comm_B, comm_union ) );
    }
    else
    {
        TEST_ASSERT( CommTools::equal( comm_A, comm_union ) );
        TEST_ASSERT( !CommTools::equal( comm_B, comm_union ) );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommTools, union_test_5 )
{
    using namespace DataTransferKit;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm comm_A = getDefaultComm<int>();
    int inverse_rank = comm_A->getSize() - comm_A->getRank() - 1;
    RCP_Comm comm_B = comm_A->split( 0, inverse_rank);
    RCP_Comm comm_global = comm_A;

    RCP_Comm comm_union;
    CommTools::unite( comm_A, comm_B, comm_union, comm_global );
    TEST_ASSERT( CommTools::equal( comm_A, comm_union ) );
    TEST_ASSERT( CommTools::equal( comm_B, comm_union ) );
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

    if ( comm_A->getSize() > 1 )
    {
        TEST_ASSERT( !CommTools::equal( comm_A, comm_intersect ) );
    }
    else
    {
        TEST_ASSERT( CommTools::equal( comm_A, comm_intersect ) );
    }

    TEST_ASSERT( CommTools::equal( comm_B, comm_intersect ) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommTools, intersect_test_3 )
{
    using namespace DataTransferKit;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm comm_A = getDefaultComm<int>();
    int inverse_rank = comm_A->getSize() - comm_A->getRank() - 1;
    RCP_Comm comm_B = comm_A->split( 0, inverse_rank);

    RCP_Comm comm_intersect;
    CommTools::intersect( comm_A, comm_B, comm_intersect );
    TEST_ASSERT( CommTools::equal( comm_A, comm_intersect ) );
    TEST_ASSERT( CommTools::equal( comm_B, comm_intersect ) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CommTools, intersect_test_4 )
{
    using namespace DataTransferKit;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > RCP_Comm;

    RCP_Comm comm_A = getDefaultComm<int>();
    int inverse_rank = comm_A->getSize() - comm_A->getRank() - 1;
    RCP_Comm comm_B = comm_A->split( 0, inverse_rank);
    RCP_Comm comm_global = comm_A;

    RCP_Comm comm_intersect;
    CommTools::intersect( comm_A, comm_B, comm_intersect, comm_global );
    TEST_ASSERT( CommTools::equal( comm_A, comm_intersect ) );
    TEST_ASSERT( CommTools::equal( comm_B, comm_intersect ) );
}

//---------------------------------------------------------------------------//
// end tstCommTools.cpp
//---------------------------------------------------------------------------//
