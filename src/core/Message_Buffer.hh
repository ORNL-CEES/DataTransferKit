//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Message_Buffer.hh
 * \author Stuart Slattery
 * \date   Tue Nov 08 07:35:21 2011
 * \brief  Message_Buffer class definition
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef core_Message_Buffer_hh
#define core_Message_Buffer_hh

#include "comm/global.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class Message_Buffer
 * \brief A simple container for buffers to be sent in MPI messages.
 */
//===========================================================================//

template<class OrdinateType_T>
class Message_Buffer 
{
  public:

    //@{
    //! Useful Typedefs.
    typedef OrdinateType_T              OrdinateType;
    typedef nemesis::Request            Request;
    typedef std::vector<char>           Buffer;
    //@}

  private:

    // MPI ordinate.
    OrdinateType d_ordinate;

    // Buffer for storing POD.
    Buffer d_buffer;

    // Request handle.
    Request d_request;

  public:
    //! Constructor.
    explicit Receive_Buffer(OrdinateType ordinate, OrdinateType buffer_size)
        : d_ordinate(ordinate)
        , d_buffer(buffer_size)
    { /* ... */ }

    //! Get the ordinate.
    OrdinateType ordinate() const { return d_ordinate; }

    //! Get the requests.
    Request& request() { return d_request; }

    //! Get the buffer.
    Buffer& buffer() { return d_buffer; }

    //! Return whether given request is complete.
    static bool complete(Receive_Buffer<OrdinateType>& buf) 
    { 
        return buf.request().complete(); 
    }
};

} // end namespace coupler

#endif // core_Message_Buffer_hh

//---------------------------------------------------------------------------//
//              end of core/Message_Buffer.hh
//---------------------------------------------------------------------------//
