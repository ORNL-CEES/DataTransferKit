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
 * \file DTK_BufferCommunicator_impl.hpp
 * \author Stuart R. Slattery
 * \brief BufferCommunicator class implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BUFFERCOMMUNICATOR_IMPL_HPP
#define DTK_BUFFERCOMMUNICATOR_IMPL_HPP

#include "DTK_DBC.hpp"

#include <Teuchos_as.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
//---------------------------------------------------------------------------//
template<class T>
BufferCommunicator<T>::BufferCommunicator( 
    const Teuchos::Array<int>& send_ranks,
    const Teuchos::Array<int>& receive_ranks,
    const Teuchos::RCP<const Comm>& comm,
    const int max_buffer_size )
    : d_send_ranks( send_ranks )
    , d_receive_ranks( receive_ranks )
    , d_num_send_neighbors( d_send_ranks.size() )
    , d_num_receive_neighbors( d_receive_ranks.size() )
    , d_size( comm->getSize() )
    , d_rank( comm->getRank() )
    , d_sends( d_num_send_neighbors )
    , d_receives( d_num_receive_neighbors )
{
    DTK_REQUIRE( Teuchos::nonnull(comm) );
    DTK_REQUIRE( d_num_send_neighbors >= 0 );
    DTK_REQUIRE( d_num_receive_neighbors >= 0 );

    // Set the static packet state.
    BDT::setByteSize();
    DTK_INSIST( BDT::getPackedBytes() );
    DataBufferType::setSizePackedData( BDT::getPackedBytes() );

    // Set the max number of packets that will be stored in each buffer.
    DataBufferType::setMaxNumPackets( max_buffer_size );

    // Duplicate the input communicator so we have a blocking and nonblocking
    // tag for the send and receive buffers.
    Teuchos::RCP<const Comm> comm_blocking = comm->duplicate();
    Teuchos::RCP<const Comm> comm_nonblocking = comm->duplicate();

    // Allocate the send buffers and set their communicators.
    for ( int n = 0; n < d_num_send_neighbors; ++n )
    {
	d_sends[n].setComm( comm_blocking, comm_nonblocking );
	d_sends[n].allocate();
    }

    // Allocate the receive buffers and set their communicators.
    for ( int n = 0; n < d_num_receive_neighbors; ++n )
    {
	d_receives[n].setComm( comm_blocking, comm_nonblocking );
	d_receives[n].allocate();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Buffer and send a packet.
 */
template<class T>
const typename BufferCommunicator<T>::Result& 
BufferCommunicator<T>::communicate( 
    const Teuchos::RCP<T>& packet, const int neighbor_id )
{
    DTK_REQUIRE( Teuchos::nonnull(packet) );

    // Initialize result status.
    d_result.sent = false;
    d_result.destination = 0;

    // Add the data to the appropriate buffer.
    d_sends[neighbor_id].bufferData( *packet );

    // Update the result destination.
    d_result.destination = d_send_ranks[neighbor_id];
    DTK_CHECK( d_result.destination < d_size );

    // If the buffer is full send it.
    if ( d_sends[neighbor_id].isFull() )
    {
	DTK_CHECK( d_sends[neighbor_id].numPackets() == 
	       Teuchos::as<int>(maxBufferSize()) );

	d_sends[neighbor_id].post( d_result.destination );
	d_sends[neighbor_id].wait();

	DTK_CHECK( d_sends[neighbor_id].isEmpty() );
	DTK_CHECK( d_sends[neighbor_id].allocatedSize() > 0 );
	DTK_CHECK( !d_sends[neighbor_id].status() );

	d_result.sent = true;
    }

    // Return the result.
    return d_result;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send all buffers that are not empty.
 */
template<class T>
int BufferCommunicator<T>::send()
{
    int num_sent = 0;

    for ( int n = 0; n < d_num_send_neighbors; ++n )
    {
	DTK_CHECK( d_sends[n].allocatedSize() > 0 );
	DTK_CHECK( d_send_ranks[n] < d_size );

	if( !d_sends[n].isEmpty() )
	{
	    DTK_CHECK( d_sends[n].numPackets() > 0 );

	    num_sent += d_sends[n].numPackets();
	    d_sends[n].post( d_send_ranks[n] );
	    d_sends[n].wait();

	    DTK_CHECK( num_sent > 0 );
	}

	DTK_ENSURE( d_sends[n].isEmpty() );
	DTK_ENSURE( d_sends[n].allocatedSize() > 0 );
	DTK_ENSURE( !d_sends[n].status() );
    }

    return num_sent;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send all buffers whether they are empty or not.
 */
template<class T>
int BufferCommunicator<T>::flush()
{
    int num_sent = 0;

    for ( int n = 0; n < d_num_send_neighbors; ++n )
    {
	DTK_CHECK( d_sends[n].allocatedSize() > 0 );
	DTK_CHECK( d_send_ranks[n] < d_size );

	num_sent += d_sends[n].numPackets();
	d_sends[n].post( d_send_ranks[n] );
	d_sends[n].wait();

	DTK_ENSURE( d_sends[n].isEmpty() );
	DTK_ENSURE( d_sends[n].allocatedSize() > 0 );
	DTK_ENSURE( !d_sends[n].status() );
    }

    return num_sent;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Post receives.
 */
template<class T>
void BufferCommunicator<T>::post()
{
    for ( int n = 0; n < d_num_receive_neighbors; ++n )
    {
	DTK_CHECK( !d_receives[n].status() );
	DTK_CHECK( d_receives[n].allocatedSize() > 0 );
	DTK_CHECK( d_receives[n].isEmpty() );
	DTK_CHECK( d_receive_ranks[n] < d_size );

	d_receives[n].post( d_receive_ranks[n] );

	DTK_ENSURE( d_receives[n].status() );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wait on receive buffers.
 */
template<class T>
int BufferCommunicator<T>::wait( BankType& bank )
{
    int num_received = 0;

    for ( int n = 0; n < d_num_receive_neighbors; ++n )
    {
	DTK_CHECK( d_receives[n].allocatedSize() > 0 );

	d_receives[n].wait();
	num_received += d_receives[n].numPackets();
	d_receives[n].addToBank( bank );

	DTK_ENSURE( !d_receives[n].status() );
	DTK_ENSURE( d_receives[n].isEmpty() );
    }

    return num_received;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Receive buffers and repost.
 */
template<class T>
int BufferCommunicator<T>::checkAndPost( BankType& bank )
{
    int num_received = 0;

    for ( int n = 0; n < d_num_receive_neighbors; ++n )
    {
	DTK_CHECK( d_receives[n].allocatedSize() > 0 );
	DTK_CHECK( d_receive_ranks[n] < d_size );

	if( d_receives[n].check() )
	{
	    num_received += d_receives[n].numPackets();
	    d_receives[n].addToBank( bank );

	    DTK_CHECK( d_receives[n].isEmpty() );
	    DTK_CHECK( !d_receives[n].status() );

	    d_receives[n].post( d_receive_ranks[n] );

	    DTK_ENSURE( d_receives[n].status() );
	}
    }

    return num_received;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Status of send buffers.
 *
 * Return true only if all send buffers are on.
 */
template<class T>
bool BufferCommunicator<T>::sendStatus()
{
    if ( d_num_send_neighbors == 0 ) return false;
      
    for ( int n = 0; n < d_num_send_neighbors; ++n )
    {
	DTK_CHECK( d_sends[n].allocatedSize() > 0 );

	if ( !d_sends[n].status() ) 
	{
	    return false;
	}
    }

    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Status of receive buffers.
 *
 * Return true only if all receive buffers are on.
 */
template<class T>
bool BufferCommunicator<T>::receiveStatus()
{
    if ( d_num_receive_neighbors == 0 ) return false;

    for ( int n = 0; n < d_num_receive_neighbors; ++n )
    {
	DTK_CHECK( d_receives[n].allocatedSize() > 0 );

	if ( !d_receives[n].status() ) 
	{
	    return false;
	}
    }

    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief End communication. Send all buffers and receive them.
 *
 * This will clear all communication requests. All buffers will be emptied and
 * therefore all data they contain lost.
 *
 * All receives must be posted or the flush will hang.
 */
template<class T>
void BufferCommunicator<T>::end()
{
    flush();

    for ( int n = 0; n < d_num_receive_neighbors; ++n )
    {
	DTK_CHECK( d_receives[n].allocatedSize() > 0 );

	d_receives[n].wait();
	d_receives[n].empty();

	DTK_ENSURE( d_receives[n].isEmpty() );
    }

    DTK_ENSURE( !sendStatus() );
    DTK_ENSURE( !receiveStatus() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Number of packets in all buffers.
 */
template<class T>
std::size_t BufferCommunicator<T>::sendBufferSize() const
{
    int send_num = 0;

    for ( int n = 0; n < d_num_send_neighbors; ++n )
    {
	DTK_CHECK( d_sends[n].allocatedSize() > 0 );

	send_num += d_sends[n].numPackets();
    }

    return send_num;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_BUFFERCOMMUNICATOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_BufferCommunicator_impl.hpp
//---------------------------------------------------------------------------//

