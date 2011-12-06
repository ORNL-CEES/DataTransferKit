//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstMessage_Buffer.cc
 * \author Stuart Slattery
 * \date   Tue Nov 08 12:20:24 2011
 * \brief  Unit test for the message buffer container.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "comm/global.hh"
#include "../Packing_Utils.hh"
#include "../Message_Buffer.hh"

using denovo::Packer;
using denovo::Unpacker;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

namespace coupler {

TEUCHOS_UNIT_TEST( Message_Buffer, messagebufferCommunication )
{
    // Make a message buffer.
    Message_Buffer<int> message_buffer(nemesis::node(), sizeof(int));

    // Check the initialization.
    TEST_ASSERT( message_buffer.ordinate() == nemesis::node() );
    TEST_ASSERT( message_buffer.buffer().size() == sizeof(int) );

    // Setup an MPI_Request for a non-blocking receive.
    nemesis::receive_async(message_buffer.request(),
			   &message_buffer.buffer()[0],
			   message_buffer.buffer().size(),
			   message_buffer.ordinate());

    // Every node sends a buffer to itself containing its node id.
    std::vector<char> send_buffer( sizeof(int) );
    Packer p;
    p.set_buffer( send_buffer.size(), &send_buffer[0] );
    p << nemesis::node();

    nemesis::send_async(&send_buffer[0], send_buffer.size(), nemesis::node());

    // Wait for the request to be processed.
    while ( !message_buffer.request().complete() ) { }

    // Unpack the buffer and check it.
    int data;
    Unpacker u;
    u.set_buffer( message_buffer.buffer().size(), 
		  &message_buffer.buffer()[0] );
    u >> data;
    
    TEST_ASSERT( data == nemesis::node() );
}

} // end namespace coupler

//---------------------------------------------------------------------------//
//                        end of tstMessage_Buffer.cc
//---------------------------------------------------------------------------//
