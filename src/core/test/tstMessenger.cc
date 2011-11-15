//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/test/tstMessenger.cc
 * \author Stuart Slattery
 * \date   Thu Jun 02 09:10:58 2011
 * \brief  Tests the Messenger class.
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
#include "utils/Packing_Utils.hh"
#include "../Messenger.hh"

using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

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

// messenger test
void messenger_test(Parallel_Unit_Test &ut)
{

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "Messenger test ok on " << nemesis::node();
        ut.passes( m.str() );
    }

}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Parallel_Unit_Test ut(argc, argv, coupler::release);

    try
    {
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;

        if(nemesis::nodes() > 1)
        {
            messenger_test(ut);
            gpass += ut.numPasses;
            gfail += ut.numFails;
            ut.reset();
        }
        else
        {
            ++gpass;
        }

        // add up global passes and fails
        nemesis::global_sum(gpass);
        nemesis::global_sum(gfail);
        ut.numPasses = gpass;
        ut.numFails  = gfail;
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstMessenger, " 
                  << err.what()
                  << std::endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstMessenger, " 
                  << "An unknown exception was thrown."
                  << std::endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstMessenger.cc
//---------------------------------------------------------------------------//
