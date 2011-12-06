//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/test/tstLG_Indexer.cc
 * \author Stuart R. Slattery
 * \date   Thu Jun 16 17:00:12 2011
 * \brief  LG_Indexer unit tests
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template_c4_test.cc,v 1.7 2008/01/02 22:50:26 9te Exp $
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "comm/global.hh"
#include "comm/Parallel_Unit_Test.hh"
#include "release/Release.hh"
#include "../LG_Indexer.hh"

#include "Teuchos_RCP.hpp"

using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

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

void indexer_test(Parallel_Unit_Test &ut)
{
    typedef coupler::LG_Indexer        LG_Indexer_t;

    // Create the test app
    Teuchos::RCP<test_app> app_ptr( new test_app() );

    // Create the local communicator object
    nemesis::Communicator_t local_comm;

    // Split the communicator
    nemesis::split(0, nodes-node-1, local_comm);

    // Make the indexer
    LG_Indexer_t indexer(get_comm_world(), local_comm, app_ptr);

    // check the map
    UNIT_TEST(indexer.size() == nodes);
    for (int i = 0; i < nodes; ++i)
    {
        UNIT_TEST(indexer.l2g(i) == nodes - i - 1);
    }

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "LG_Indexer test passes on " << node;
        ut.passes( m.str() );
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Parallel_Unit_Test ut(argc, argv, coupler::release);

    node  = nemesis::node();
    nodes = nemesis::nodes();
    
    try
    {
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;
 
        indexer_test(ut);
        gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();
       
        // add up global passes and fails
        nemesis::global_sum(gpass);
        nemesis::global_sum(gfail);
        ut.numPasses = gpass;
        ut.numFails  = gfail;
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstLG_Indexer, " 
                  << err.what()
                  << std::endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstLG_Indexer, " 
                  << "An unknown exception was thrown."
                  << std::endl;
        ut.numFails++;
    }
    return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstLG_Indexer.cc
//---------------------------------------------------------------------------//
