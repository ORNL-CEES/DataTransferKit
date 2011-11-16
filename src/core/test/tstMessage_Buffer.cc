//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstMessage_Buffer.cc
 * \author Stuart Slattery
 * \date   Tue Nov 08 12:20:24 2011
 * \brief  Unit test for the message buffer container.
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
#include "../Message_Buffer.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

using denovo::Packer;
using denovo::Unpacker;

using coupler::Message_Buffer;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void message_buffer_test(Parallel_Unit_Test &ut)
{
    // Make a message buffer.
    Message_Buffer<int> message_buffer(node, sizeof(int));

    // Check the initialization.
    UNIT_TEST( message_buffer.ordinate() == node );
    UNIT_TEST( message_buffer.buffer().size() == sizeof(int) );

    // Setup an MPI_Request for a non-blocking receive.
    nemesis::receive_async(message_buffer.request(),
			   &message_buffer.buffer()[0],
			   message_buffer.buffer().size(),
			   message_buffer.ordinate());

    // Every node sends a buffer to itself containing its node id.
    std::vector<char> send_buffer( sizeof(int) );
    Packer p;
    p.set_buffer( send_buffer.size(), &send_buffer[0] );
    p << node;

    nemesis::send_async(&send_buffer[0], send_buffer.size(), node);

    // Wait for the request to be processed.
    while ( !message_buffer.request().complete() ) { }

    // Unpack the buffer and check it.
    int data;
    Unpacker u;
    u.set_buffer( message_buffer.buffer().size(), 
		  &message_buffer.buffer()[0] );
    u >> data;
    
    UNIT_TEST( data == node );

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "Transfer_Map test passes on " << node;
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
        
	message_buffer_test(ut);
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
        std::cout << "ERROR: While testing tstMessage_Buffer, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstMessage_Buffer, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstMessage_Buffer.cc
//---------------------------------------------------------------------------//
