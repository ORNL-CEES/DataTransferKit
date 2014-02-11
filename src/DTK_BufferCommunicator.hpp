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
 * \file DTK_BufferCommunicator.hpp
 * \author Stuart R. Slattery
 * \brief BufferCommunicator class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BUFFERCOMMUNICATOR_HPP
#define DTK_BUFFERCOMMUNICATOR_HPP

#include "DTK_BufferDataTraits.hpp"
#include "DTK_DataBuffer.hpp"
#include "DTK_CommDataBuffer.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class BufferCommunicator 
 * \brief Structure for communicating data packets.
 *
 * Tom Evans is responsible for the design of this class.
 */
//---------------------------------------------------------------------------//
template<class T>
class BufferCommunicator
{
  public:

    //@{
    //! Typedefs.
    typedef T                                            data_type;
    typedef BufferDataTraits<T>                          BDT;
    typedef DataBuffer<T>                                DataBufferType;
    typedef typename DataBufferType::BankType            BankType;
    typedef SendDataBuffer<T>                            SendBuffer;
    typedef ReceiveDataBuffer<T>                         ReceiveBuffer;
    typedef Teuchos::Comm<int>                           Comm;
    //@}

    //! Communication result.
    struct Result
    {
        bool sent;
        int  destination;
    };

  public:

    // Constructor.
    BufferCommunicator( const Teuchos::Array<int>& send_ranks,
			const Teuchos::Array<int>& receive_ranks,
			const Teuchos::RCP<const Comm>& comm, 
			const int max_buffer_size );

    // Destructor.
    ~BufferCommunicator()
    { /* ... */ }

    // Buffer and send a packet.
    const Result& communicate( const Teuchos::RCP<T>& packet,
			       const int neighbor_id );

    // Send all buffers that are not empty.
    int send();

    // Flush all buffers whether they are empty or not.
    int flush();

    // Post receives.
    void post();

    // Wait on receive buffers.
    int wait( BankType& bank );

    // Receive buffers and repost.
    int checkAndPost( BankType& bank );

    // Status of send buffers.
    bool sendStatus();

    // Status of receive buffers.
    bool receiveStatus();

    // End communication.
    void end();

    //! Packet buffer size.
    std::size_t maxBufferSize() const
    { return DataBufferType::maxNum(); }

    // Number of packets in all buffers.
    std::size_t sendBufferSize() const;

    // Get a send buffer by local id.
    const SendBuffer& sendBuffer( int n ) const
    { return d_sends[n]; }

    // Get a receive buffer by local id.
    const ReceiveBuffer& receiveBuffer( int n ) const
    { return d_receives[n]; }

  private:

    // The neighbors we are sending to.
    Teuchos::Array<int> d_send_ranks;

    // The neighbors we are receiving from.
    Teuchos::Array<int> d_receive_ranks;

    // Number of neighbors we are sending to.
    int d_num_send_neighbors;

    // Number of neighbors we are receiving from.
    int d_num_receive_neighbors;
    // Communicator size.
    int d_size;

    // Communicator rank.
    int d_rank;

    // Send buffers.
    Teuchos::Array<SendBuffer> d_sends;

    // Receive buffers.
    Teuchos::Array<ReceiveBuffer> d_receives;

    // Result of a packet communication.
    Result d_result;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_BufferCommunicator_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_BUFFERCOMMUNICATOR_HPP

//---------------------------------------------------------------------------//
// end DTK_BufferCommunicator.hpp
//---------------------------------------------------------------------------//

