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
 * \file DTK_DataBuffer_impl.hpp
 * \author Stuart R. Slattery
 * \brief DataBuffer class implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATABUFFER_IMPL_HPP
#define DTK_DATABUFFER_IMPL_HPP

#include <algorithm>

#include "DTK_DBC.hpp"
#include "DTK_Serializer.hpp"

#include <Teuchos_as.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Size constructor.
 */
template<class T>
DataBuffer<T>::DataBuffer( std::size_t size, int num_packes )
    : d_number( 0 )
{
    setSizePackedData( size );
    setMaxNumPackets( num_packets );
    allocate();
    DTK_ENSURE( !d_buffer.empty() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Allocate the buffer for the byte size of the maximum number of
 * packets plus an additional integer for the actual number of the buffer.
 */
template<class T>
void DataBuffer<T>::allocate()
{
    DTK_REQUIRE( d_number == 0 );
    d_buffer.resize( 
	d_max_num_packets*d_size_packed_data + sizeof(int), '\0' );
    DTK_ENSURE( d_number == 0 );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Deallocate the buffer.
 */
template<class T>
void DataBuffer<T>::deallocate()
{
    DTK_REQUIRE( d_number == 0 );
    d_buffer.clear();
    DTK_ENSURE( d_number == 0 );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write a packet into the buffer.
 */
template<class T>
void DataBuffer<T>::bufferData( const T& packet )
{
    DTK_REQUIRE( d_size_packed_data > 0 );
    DTK_REQUIRE( d_number < d_max_num_packets );
    DTK_REQUIRE( d_number >= 0 );
    DTK_REQUIRE( !d_buffer.empty() );

    Buffer packed_data = BDT::pack( data );
    DTK_CHECK( Teuchos::as<std::size_t>(packed_data.size()) == 
	       d_size_packed_data );

    Buffer::iterator buffer_it = d_buffer.begin() + 
				 d_size_packed_data*d_number;
    DTK_REQUIRE( buffer_it != d_buffer.end() );

    std::copy( packed_data.begin(), packed_data.end(), buffer_it );
    ++d_number;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add the packets in the buffer to a bank.
 */
template<class T>
void DataBuffer<T>::addToBank( BankType& bank )
{
    DTK_REQUIRE( d_size_packed_data > 0 );

    Buffer::const_iterator buffer_it = d_buffer.begin();
    Buffer packed_data( d_size_packed_data );
    Teuchos::RCP<T> data;

    DTK_REMEMBER( std::size_t bank_size = bank.size() );

    for ( int n = 0; n < d_number; ++n )
    {
	std::copy( buffer_it, buffer_it + d_size_packed_data, 
		   packed_data.begin() );

	data = BDT::createFromBuffer( packed_data );
	DTK_CHECK( !data.is_null() );
	bank.push( data );

	buffer_it += d_size_packed_data;
    }

    DTK_ENSURE( bank_size + d_number == bank.size() );
    DTK_ENSURE( d_number == d_max_num_packets ?
		buffer_it + sizeof(int) == d_buffer.end() :
		buffer_it + sizeof(int) != d_buffer.end() );

    empty();
    DTK_ENSURE( isEmpty() );
}

//---------------------------------------------------------------------------//
// Protected Members.
//---------------------------------------------------------------------------//
/*!
 * \brief Add the number of packets to the end of the buffer.
 */
template<class T>
void DataBuffer<T>::writeNumToBuffer()
{
    DTK_REQUIRE( Teuchos::as<std::size_t>(d_buffer.size()) > sizeof(int) );
    Serializer s;
    s.setBuffer( sizeof(int), &d_buffer[d_buffer.size() - sizeof(int)] );
    s << d_number;
    DTK_ENSURE( s.getPtr() == &d_buffer[0] + d_buffer.size() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read the number of packets from the end of the buffer.
 */
template<class T>
void DataBuffer<T>::readNumFromBuffer()
{
    Deserializer ds;
    ds.setBuffer( sizeof(int), &d_buffer[d_buffer.size() - sizeof(int)] );
    ds >> d_number;
    DTK_ENSURE( ds.getPtr() == &d_buffer[0] + d_buffer.size() );
    DTK_ENSURE( d_number >= 0 );
}

//---------------------------------------------------------------------------//
// Static Members.
//---------------------------------------------------------------------------//

//! Default maximum number of data packets allowed in a buffer.
template<class T>
int DataBuffer<T>::d_max_num_packets = 1000;

//! Default size of a packed data packet.
template<class T>
std::size_t DataBuffer<T>::d_size_packed_data = 0;

//---------------------------------------------------------------------------//
/*!
 * \brief Set the maximum number of data packets allowed in the buffer.
 */
template<class T>
void DataBuffer<T>::setMaxNumPackets( int num_packets )
{
    DTK_REQUIRE( num_packets > 0 );
    DTK_REQUIRE( d_size_packed_data > 0 );
    d_max_num_packets = num_packets;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set the byte size of a packed data packet.
 */
template<class T>
void DataBuffer<T>::setSizePackedData( std::size_t size )
{
    DTK_REQUIRE( size > 0 );
    d_size_packed_data = size;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_DATABUFFER_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_DataBuffer_impl.hpp
//---------------------------------------------------------------------------//

