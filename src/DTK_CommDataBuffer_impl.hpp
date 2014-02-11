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
 * \file DTK_CommDataBuffer_impl.hpp
 * \author Stuart R. Slattery
 * \brief CommDataBuffer class implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_COMMDATABUFFER_IMPL_HPP
#define DTK_COMMDATABUFFER_IMPL_HPP

#include <DTK_CommTools.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Ptr.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// CommDataBuffer functions.
//---------------------------------------------------------------------------//
/*!
 * \brief Pure virtual destructor. Prevents direct instantiation of
 * CommDataBuffer.
 */
template<class T>
CommDataBuffer<T>::~CommDataBuffer()
{ /* ... */ }

//---------------------------------------------------------------------------//
// ReceiveDataBuffer functions.
//---------------------------------------------------------------------------//
/*!
 * \brief Blocking receive.
 */
template<class T>
void ReceiveDataBuffer<T>::receive( int rank )
{
    DTK_REQUIRE( !Base::d_comm_blocking.is_null() );
    DTK_REQUIRE( Root::isEmpty() );
    DTK_REQUIRE( Root::allocatedSize() > sizeof(int) );

    Teuchos::receive<int,char>( 
	*Base::d_comm_blocking, rank, 
	Root::d_buffer.size(), Root::d_buffer.getRawPtr() );
    Root::readNumFromBuffer();

    DTK_ENSURE( Root::d_number < Root::maxNum() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Post non-blocking receives.
 */
template<class T>
void ReceiveDataBuffer<T>::post( int rank )
{
    DTK_REQUIRE( !Base::d_comm_nonblocking.is_null() );
    DTK_REQUIRE( Root::isEmpty() );
    DTK_REQUIRE( Root::allocatedSize() > sizeof(int) );

    Base::d_handle = Teuchos::ireceive<int,char>( 
	*Base::d_comm_nonblocking, 
	Teuchos::arcpFromArray(Root::d_buffer), 
	rank );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wait on a non-blocking receive to finish.
 */
template<class T>
void ReceiveDataBuffer<T>::wait()
{
    DTK_REQUIRE( !Base::d_comm_nonblocking.is_null() );

    Teuchos::Ptr<Teuchos::RCP<typename Base::Request> > 
	request_ptr( &this->d_handle );
    Teuchos::wait( *Base::d_comm_nonblocking, request_ptr );
    Root::readNumFromBuffer();

    DTK_ENSURE( Base::d_handle.is_null() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check to see if a non-blocking send has finished.
 */
template<class T>
bool ReceiveDataBuffer<T>::check()
{
    DTK_REQUIRE( !Base::d_comm_nonblocking.is_null() );

    if ( CommTools::isRequestComplete(Base::d_handle) )
    {
	Root::readNumFromBuffer();
	Base::d_handle = Teuchos::null;

	DTK_ENSURE( Base::d_handle.is_null() );
	DTK_ENSURE( Root::numHistories() >= 0 );
	return true;
    }

    return false;
}

//---------------------------------------------------------------------------//
// SendDataBuffer functions.
//---------------------------------------------------------------------------//
/*!
 * \brief Blocking send.
 */
template<class T>
void SendDataBuffer<T>::send( int rank )
{
    DTK_REQUIRE( !Base::d_comm_blocking.is_null() );
    DTK_REQUIRE( Root::allocatedSize() > sizeof(int) );

    Root::writeNumToBuffer();
    Teuchos::send<int,char>( *Base::d_comm_blocking, Root::d_buffer.size(), 
			     Root::d_buffer.getRawPtr(), rank );

    Root::empty();

    DTK_ENSURE( Root::isEmpty() );
    DTK_ENSURE( Root::allocatedSize() > sizeof(int) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Post non-blocking send.
 */
template<class T>
void SendDataBuffer<T>::post( int rank )
{
    DTK_REQUIRE( !Base::d_comm_nonblocking.is_null() );
    DTK_REQUIRE( Root::allocatedSize() > sizeof(int) );

    Root::writeNumToBuffer();
    Base::d_handle = Teuchos::isend<int,char>( 
	*Base::d_comm_nonblocking, 
	Teuchos::arcpFromArray(Root::d_buffer), 
	rank );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wait on a non-blocking send to finish.
 */
template<class T>
void SendDataBuffer<T>::wait()
{
    DTK_REQUIRE( !Base::d_comm_nonblocking.is_null() );

    Teuchos::Ptr<Teuchos::RCP<typename Base::Request> > 
	request_ptr( &this->d_handle );

    Teuchos::wait( *Base::d_comm_nonblocking, request_ptr );

    Root::empty();

    DTK_ENSURE( Base::d_handle.is_null() );
    DTK_ENSURE( Root::isEmpty() );
    DTK_ENSURE( Root::allocatedSize() > sizeof(int) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check to see if a non-blocking send has finished.
 */
template<class T>
bool SendDataBuffer<T>::check()
{
    DTK_REQUIRE( !Base::d_comm_nonblocking.is_null() );

    if ( CommTools::isRequestComplete(Base::d_handle) )
    {
	Root::empty();
	Base::d_handle = Teuchos::null;

	DTK_ENSURE( Base::d_handle.is_null() );
	DTK_ENSURE( Root::isEmpty() );
	return true;
    }

    return false;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_COMMDATABUFFER_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_CommDataBuffer_impl.hpp
//---------------------------------------------------------------------------//

