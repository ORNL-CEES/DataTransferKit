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

  *: Neither the name of the Oak Ridge National Laboratory names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.

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
 * \file DTK_DataSerializer.hpp
 * \author Stuart R. Slattery
 * \brief Plain old data serializer.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DATASERIALIZER_HPP
#define DTK_DATASERIALIZER_HPP

#include <cstring>

#include "DTK_DBC.hpp"

#include <Teuchos_ArrayView.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class DataSerializer
 * \brief DataSerializer for putting data into a byte stream. Based on Tom Evan's
 * Packer class.
 */
class DataSerializer
{
  public:

    //@{
    //! Typedefs.
    typedef char               data_type;
    typedef char*              ptr_type;
    typedef const char*        const_ptr_type;
    //@}

  public:

    //! Constructor.
    DataSerializer()
	: d_size( 0 )
	, d_ptr( 0 )
	, d_begin( 0 )
	, d_end( 0 )
	, d_size_mode( false )
    { /* ... */ }

    //! Destructor.
    ~DataSerializer()
    { /* ... */ }

    //! Set the buffer and put into pack mode.
    inline void setBuffer( const std::size_t& size, ptr_type buffer )
    {
	DTK_REQUIRE( buffer );
	d_size_mode = false;
	d_size = size;
	d_ptr = buffer;
	d_begin = buffer;
	d_end = d_begin + d_size;
    }

    //! Set the buffer and put into pack mode.
    inline void setBuffer( const Teuchos::ArrayView<data_type>& buffer )
    {
	DTK_REQUIRE( buffer.getRawPtr() );
	d_size_mode = false;
	d_size = buffer.size();
	d_ptr = buffer.getRawPtr();
	d_begin = buffer.getRawPtr();
	d_end = d_begin + d_size;
    }

    //! Put into compute buffer size mode.
    void computeBufferSizeMode()
    {
	d_size = 0;
	d_size_mode = true;
    }

    //! Pack values into the buffer.
    template<class T>
    inline void pack( const T& data )
    {
	if ( d_size_mode )
	    d_size += sizeof(T);
	else
	{
	    DTK_REQUIRE( d_begin );
	    DTK_REQUIRE( d_ptr >= d_begin);
	    DTK_REQUIRE( d_ptr + sizeof(T) <= d_end );
	
	    std::memcpy( d_ptr, &data, sizeof(T) );
	
	    d_ptr += sizeof(T);
	}
    }

    //! Get a pointer to the current position of the data stream.
    const_ptr_type getPtr() const 
    { 
	DTK_REQUIRE( !d_size_mode ); 
	return d_ptr; 
    }

    //! Get a pointer to the beginning position of the data stream.
    const_ptr_type begin() const 
    { 
	DTK_REQUIRE( !d_size_mode ); 
	return d_begin; 
    }

    //! Get a pointer to the ending position of the data stream.
    const_ptr_type end() const 
    { 
	DTK_REQUIRE( !d_size_mode ); 
	return d_end; 
    }

    //! Get the size of the data stream.
    std::size_t size() const 
    { 
	return d_size; 
    }
    
  private:

    // Size of buffer.
    std::size_t d_size;

    // Pointer to current location in buffer.
    ptr_type d_ptr;

    // Pointer the the beginning of the buffer.
    ptr_type d_begin;

    // Pointer to the end of the buffer.
    ptr_type d_end;

    // Boolean for size mode.
    bool d_size_mode;
};
    
//---------------------------------------------------------------------------//
/*!
 * \class DataDeserializer
 * \brief DataDeserializer for pulling data out of a byte stream. Based on Tom
 * Evan's Unpacker class. 
 */
class DataDeserializer
{
  public:

    //@{
    //! Typedefs.				
    typedef char               data_type;
    typedef char*              ptr_type;
    typedef const char*        const_ptr_type;
    //@}

  public:

    //! Constructor.
    DataDeserializer()
	: d_size( 0 )
	, d_ptr( 0 )
	, d_begin( 0 )
	, d_end( 0 )
    { /* ... */ }

    //! Destructor.
    ~DataDeserializer()
    { /* ... */ }

    //! Set the buffer.
    inline void setBuffer( const std::size_t& size, ptr_type buffer )
    {
	DTK_REQUIRE( buffer );
	d_size = size;
	d_ptr = buffer;
	d_begin = buffer;
	d_end = d_begin + d_size;
    }

    //! Set the buffer.
    inline void setBuffer( const Teuchos::ArrayView<data_type>& buffer )
    {
	DTK_REQUIRE( buffer.getRawPtr() );
	d_size = buffer.size();
	d_ptr = buffer.getRawPtr();
	d_begin = buffer.getRawPtr();
	d_end = d_begin + d_size;
    }

    //! Pack values into the buffer.
    template<class T>
    inline void unpack( T& data )
    {
	DTK_REQUIRE( d_begin );
	DTK_REQUIRE( d_ptr >= d_begin);
	DTK_REQUIRE( d_ptr + sizeof(T) <= d_end );
	
	std::memcpy( &data, d_ptr, sizeof(T) );
	
	d_ptr += sizeof(T);
    }

    //! Get a pointer to the current position of the data stream.
    const_ptr_type getPtr() const 
    { 
	return d_ptr; 
    }

    //! Get a pointer to the beginning position of the data stream.
    const_ptr_type begin() const 
    { 
	return d_begin; 
    }

    //! Get a pointer to the ending position of the data stream.
    const_ptr_type end() const 
    { 
	return d_end; 
    }

    //! Get the size of the data stream.
    std::size_t size() const 
    { 
	return d_size; 
    }

  private:

    // Size of buffer.
    std::size_t d_size;

    // Pointer to current location in buffer.
    ptr_type d_ptr;

    // Pointer the the beginning of the buffer.
    ptr_type d_begin;

    // Pointer to the end of the buffer.
    ptr_type d_end;
};

//---------------------------------------------------------------------------//
// Stream operators.
//---------------------------------------------------------------------------//
/*!
 * \brief Stream out (<<) operator for packing data.
 */
template<class T>
inline DataSerializer& operator<<( DataSerializer &serializer, const T &data )
{
    serializer.pack( data );
    return serializer;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Stream in (>>) operator for unpacking data.
 */
template<class T>
inline DataDeserializer& operator>>( DataDeserializer &deserializer, T &data )
{
    deserializer.unpack( data );
    return deserializer;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_DATASERIALIZER_HPP

//---------------------------------------------------------------------------//
// end DTK_DataSerializer.hpp
// ---------------------------------------------------------------------------//

