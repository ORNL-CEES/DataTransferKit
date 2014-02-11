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
 * \file DTK_CommDataBuffer.hpp
 * \author Stuart R. Slattery
 * \brief CommDataBuffer class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_COMMDATABUFFER_HPP
#define DTK_COMMDATABUFFER_HPP

#include "DTK_DBC.hpp"
#include "DTK_DataBuffer.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class CommDataBuffer
 * \brief Data buffer for communicating packets. Tom Evans is responsible
 * for the design of this class and subsequent inheritance structure.
 */
//---------------------------------------------------------------------------//
template<class T>
class CommDataBuffer : public DataBuffer<T>
{
  public:

    //@{
    //! Typedefs.
    typedef DataBuffer<T>                            Base;
    typedef typename Base::data_type                 data_type;
    typedef Teuchos::CommRequest<int>                Request;
    typedef Teuchos::Comm<int>                       Comm;
    //@}

  public:

    //! Default constructor.
    CommDataBuffer()
	: d_handle( Teuchos::null )
    { DTK_ENSURE( Base::isEmpty() ); }

    //! Comm constructor.
    CommDataBuffer( const Teuchos::RCP<const Comm>& comm_blocking,
		    const Teuchos::RCP<const Comm>& comm_nonblocking )
	: d_handle( Teuchos::null )
	, d_comm_blocking( comm_blocking )
	, d_comm_nonblocking( comm_nonblocking )
    { DTK_ENSURE( Base::isEmpty() ); }

    //! Size constructor.
    CommDataBuffer( const Teuchos::RCP<const Comm>& comm_blocking,
		    const Teuchos::RCP<const Comm>& comm_nonblocking,
		    std::size_t size, int num_packets )
	: Base( size, num_packets )
	, d_handle( Teuchos::null )
	, d_comm_blocking( comm_blocking )
	, d_comm_nonblocking( comm_nonblocking )
    {
	DTK_ENSURE( Base::isEmpty() );
	DTK_ENSURE( Base::allocatedSize() > 0 );
    }

    // Pure virtual destructor.
    virtual ~CommDataBuffer() = 0;

    //! Asynchronous post.
    virtual void post( int rank ) = 0;

    //! Asynchronous wait.
    virtual void wait() = 0;

    //! Asynchronous check.
    virtual bool check() = 0;

    //! Free non-blocking communication buffer handles.
    inline void free()
    {
	d_handle = Teuchos::null;
	Base::empty();
	DTK_ENSURE( Base::isEmpty() );
	DTK_ENSURE( d_handle.is_null() );
    }

    //! Check the status of a non-blocking communication buffer.
    inline bool status() const
    { return !d_handle.is_null(); }

    //! Set the communicators for this buffer.
    void setComm( const Teuchos::RCP<const Comm>& comm_blocking,
		  const Teuchos::RCP<const Comm>& comm_nonblocking )
    { 
	d_comm_blocking = comm_blocking; 
	d_comm_nonblocking = comm_nonblocking; 
	DTK_ENSURE( !d_comm_blocking.is_null() );
	DTK_ENSURE( !d_comm_nonblocking.is_null() );
    }

  protected:

    // Non-blocking communication handles. This object's destructor will
    // cancel the request. A handle is in use if it is non-null.
    Teuchos::RCP<Request> d_handle;

    // Communicator on which blocking communications for this buffer is
    // defined.
    Teuchos::RCP<const Comm> d_comm_blocking;

    // Communicator on which nonblocking communications for this buffer is
    // defined.
    Teuchos::RCP<const Comm> d_comm_nonblocking;
};

//---------------------------------------------------------------------------//
/*!
 * \class ReceiveDataBuffer
 * \brief Data buffer for receiving packets. Tom Evans is responsible for
 * the design of this class and subsequent inheritance structure.
 */
//---------------------------------------------------------------------------//
template<class T>
class ReceiveDataBuffer : public CommDataBuffer<T>
{
  public:

    //@{
    //! Typedefs.
    typedef DataBuffer<T>                      Root;
    typedef CommDataBuffer<T>                  Base;
    typedef typename Base::Comm                Comm;
    //@}

  public:

    //! Default constructor.
    ReceiveDataBuffer()
    { DTK_ENSURE( Base::isEmpty() ); }

    //! Comm constructor.
    ReceiveDataBuffer( const Teuchos::RCP<const Comm>& comm_blocking,
		       const Teuchos::RCP<const Comm>& comm_nonblocking )
	: Base( comm_blocking, comm_nonblocking )
    { DTK_ENSURE( Base::isEmpty() ); }

    //! Size constructor.
    ReceiveDataBuffer( const Teuchos::RCP<const Comm>& comm_blocking,
		       const Teuchos::RCP<const Comm>& comm_nonblocking,
		       std::size_t size, int num_packets )
	: Base( comm_blocking, comm_nonblocking, size, num_packets )
    {
	DTK_ENSURE( Base::isEmpty() );
	DTK_ENSURE( Base::allocatedSize() > 0 );
    }

    //! Destructor.
    ~ReceiveDataBuffer()
    { /* ... */ }

    // Blocking receive.
    void receive( int rank );

    // Asynchronous post.
    void post( int rank );

    // Asynchronous wait.
    void wait();

    // Asynchronous check.
    bool check();
};

//---------------------------------------------------------------------------//
/*!
 * \class SendDataBuffer
 * \brief Data buffer for sending packets. Tom Evans is responsible for the
 * design of this class and subsequent inheritance structure.
 */
//---------------------------------------------------------------------------//
template<class T>
class SendDataBuffer : public CommDataBuffer<T>
{
  public:

    //@{
    //! Typedefs.
    typedef DataBuffer<T>                      Root;
    typedef CommDataBuffer<T>                  Base;
    typedef typename Base::Comm                Comm;
    //@}

  public:

    //! Default constructor.
    SendDataBuffer()
    { DTK_ENSURE( Base::isEmpty() ); }

    //! Comm constructor.
    SendDataBuffer( const Teuchos::RCP<const Comm>& comm_blocking,
		    const Teuchos::RCP<const Comm>& comm_nonblocking )
	: Base( comm_blocking, comm_nonblocking )
    { DTK_ENSURE( Base::isEmpty() ); }

    //! Size constructor.
    SendDataBuffer( const Teuchos::RCP<const Comm>& comm_blocking,
		    const Teuchos::RCP<const Comm>& comm_nonblocking,
		    std::size_t size, int num_packets )
	: Base( comm_blocking, comm_nonblocking, size, num_packets )
    {
	DTK_ENSURE( Base::isEmpty() );
	DTK_ENSURE( Base::allocatedSize() > 0 );
    }

    //! Destructor.
    ~SendDataBuffer()
    { /* ... */ }

    // Blocking send.
    void send( int rank );

    // Asynchronous post.
    void post( int rank );

    // Asynchronous wait.
    void wait();

    // Asynchronous check.
    bool check();
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_CommDataBuffer_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_COMMDATABUFFER_HPP

//---------------------------------------------------------------------------//
// end DTK_CommDataBuffer.hpp
//---------------------------------------------------------------------------//

