//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \brief DTK_AbstractSerializer_impl.hpp
 * \author Stuart R. Slattery
 * \brief Serializable object policy interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ABSTRACTSERIALIZER_IMPL_HPP
#define DTK_ABSTRACTSERIALIZER_IMPL_HPP

#include <string>

#include "DTK_DBC.hpp"
#include "DTK_AbstractBuilder.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class AbstractSerializer
  \brief Serializer for abstract objects.
*/
//---------------------------------------------------------------------------//
//! This is an indirect serializer.
template<class Ordinal, class T>
const bool AbstractSerializer<Ordinal,T>::supportsDirectSerialization = false;

//---------------------------------------------------------------------------//
// Return the number of bytes for count objects.
template<class Ordinal, class T>
Ordinal AbstractSerializer<Ordinal,T>::fromCountToIndirectBytes( 
    const Ordinal count, const Packet buffer[])
{
    Ordinal byte_size = 0;
    for ( Ordinal i = 0; i < count; ++i )
    {
	// Bytes for the integral key.
	byte_size += sizeof(int);
	    
	// Bytes for the size of the object.
	byte_size += sizeof(std::size_t);

	// Bytes for the object.
	byte_size += ASOP::byteSize( buffer[i] );
    }
    return byte_size;
}

//---------------------------------------------------------------------------//
// Serialize to an indirect char buffer.
template<class Ordinal, class T>
void AbstractSerializer<Ordinal,T>::serialize( const Ordinal count, 
					       const Packet buffer[], 
					       const Ordinal bytes, 
					       char charBuffer[] )
{
    // Get the builder for the objects.
    Teuchos::RCP<AbstractBuilder<T> > builder = ASOP::getBuilder();

    // Serialize the objects.
    DTK_REQUIRE( fromCountToIndirectBytes(count,buffer) == bytes );
    char* buffer_pos = &charBuffer[0];
    std::size_t object_size = 0;
    int integral_key = 0;
    for ( Ordinal i = 0; i < count; ++i )
    {
	// Serialize the integral key of the object.
	integral_key = builder->getIntegralKey( ASOP::objectType(buffer[i]) );
	std::memcpy( buffer_pos, &integral_key, sizeof(int) );
	buffer_pos += sizeof(int);

	// Serialize the size of the object.
	object_size = ASOP::byteSize( buffer[i] );
	std::memcpy( buffer_pos, &object_size, sizeof(std::size_t) );
	buffer_pos += sizeof(std::size_t);
	    
	// Serialize the object.
	Teuchos::ArrayView<char> buffer_view( buffer_pos, object_size );
	ASOP::serialize( buffer[i], buffer_view );

	// Move the front of the buffer forward.
	buffer_pos += object_size;
    }
    DTK_ENSURE( &charBuffer[0] + bytes == buffer_pos );
}

//---------------------------------------------------------------------------//
// Return the number of objects for bytes of storage.
template<class Ordinal, class T>
Ordinal AbstractSerializer<Ordinal,T>::fromIndirectBytesToCount( 
    const Ordinal bytes, const char charBuffer[] )
{
    char* buffer_pos = const_cast<char*>(&charBuffer[0]);
    char* buffer_end = const_cast<char*>(&charBuffer[0] + bytes);
    Ordinal count = 0;
    std::size_t object_size = 0;
    while ( buffer_pos < buffer_end )
    {
	buffer_pos += sizeof(int);
	std::memcpy( &object_size, buffer_pos, sizeof(std::size_t) );
	buffer_pos += sizeof(std::size_t) + object_size;
	++count;
    }
    DTK_ENSURE( buffer_end == buffer_pos );
    return count;
}

//---------------------------------------------------------------------------//
// Deserialize from an indirect char buffer.
template<class Ordinal, class T>
void AbstractSerializer<Ordinal,T>::deserialize( const Ordinal bytes, 
						 const char charBuffer[], 
						 const Ordinal count, 
						 Packet buffer[] )
{
    // Get the builder for the objects.
    Teuchos::RCP<AbstractBuilder<T> > builder = ASOP::getBuilder();

    // Deserialize the objects.
    DTK_REQUIRE( fromIndirectBytesToCount(bytes,charBuffer) == count );
    char* buffer_pos = const_cast<char*>(&charBuffer[0]);
    std::size_t object_size = 0;
    int integral_key = 0;
    for ( Ordinal i = 0; i < count; ++i )
    {
	// Get the integral key for the object.
	std::memcpy( &integral_key, buffer_pos, sizeof(int) );
	buffer_pos += sizeof(int);

	// Create an object of the correct derived class.
	buffer[i] = builder->create( integral_key );

	// Get the size of the object.
	std::memcpy( &object_size, buffer_pos, sizeof(std::size_t) );
	buffer_pos += sizeof(std::size_t);

	// Deserialize the object.
	Teuchos::ArrayView<char> buffer_view( buffer_pos, object_size );
	ASOP::deserialize( buffer[i], buffer_view );

	// Move the front of the buffer forward.
	buffer_pos += object_size;
    }
    DTK_ENSURE( &charBuffer[0] + bytes == buffer_pos );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ABSTRACTSERIALIZER_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_AbstractSerializer_impl.hpp
//---------------------------------------------------------------------------//
