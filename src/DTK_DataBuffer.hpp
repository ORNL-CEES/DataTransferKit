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
 * \file DTK_DataBuffer.hpp
 * \author Stuart R. Slattery
 * \brief DataBuffer class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATABUFFER_HPP
#define DTK_DATABUFFER_HPP

#include <stack>

#include "DTK_BufferDataTraits.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class DataBuffer
 * \brief Data buffer. Tom Evans is responsible for the design
 * of this class and subsequent inheritance structure.
 */
//---------------------------------------------------------------------------//
template<class T>
class DataBuffer
{
  public:

    //@{
    //! Typedefs.
    typedef T                                  data_type;
    typedef BufferDataTraits<T>                BDT;
    typedef std::stack<Teuchos::RCP<T> >       BankType;
    typedef Teuchos::Array<char>               Buffer;
    //@}

    //! Default constructor.
    DataBuffer()
	: d_number( 0 )
    { /* ... */ }

    // Size constructor.
    DataBuffer( std::size_t size, int num_packets );

    //! Destructor.
    virtual ~DataBuffer()
    { /* ... */ }

    //! Set the number of packets in the buffer to zero.
    void empty()
    { d_number = 0; }

    // Allocate the buffer.
    void allocate();

    // Deallocate the buffer.
    void deallocate();

    // Write a packet into the buffer.
    void bufferData( const T& packet );

    // Add the packets in the buffer to a bank.
    void addToBank( BankType& bank );

    //! Get current number of packets in the buffer.
    int numPackets() const
    { return d_number; }

    //! Check if the buffer is empty.
    bool isEmpty() const
    { return ( d_number == 0 ); }

    //! Check if the buffer is full.
    bool isFull() const 
    { return ( d_number == d_max_num_packets ); }

    //! Get the current allocated size of the buffer.
    std::size_t allocatedSize() const
    { return d_buffer.size(); }

  public:

    // Set the maximum number of packets allowed in the buffer.
    static void setMaxNumPackets( int num_packets );

    // Set the byte size of a packed data packet.
    static void setSizePackedData( std::size_t size );

    //! Get the maximum number of packets allowed in the buffer.
    static int maxNum()
    { return d_max_num_packets; }

    //! Get the size of a packed data packet.
    static int sizePackedData()
    { return d_size_packed_data; }

  protected:

    // Add the number of packets to the end of the buffer.
    void writeNumToBuffer();

    // Read the number of packets from the end of the buffer.
    void readNumFromBuffer();

  protected:

    // Packed data buffer.
    Buffer d_buffer;

    // Number of packets currently in the buffer.
    int d_number;

  private:

    // Maximum number of packets allowed in the buffer.
    static int d_max_num_packets;

    // Size of a packed data in bytes.
    static std::size_t d_size_packed_data; 
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_DataBuffer_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_DATABUFFER_HPP

//---------------------------------------------------------------------------//
// end DTK_DataBuffer.hpp
//---------------------------------------------------------------------------//

