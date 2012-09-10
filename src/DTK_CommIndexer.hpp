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
 * \file DTK_CommIndexer.hpp
 * \author Stuart Slattery
 * \brief CommIndexer declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_COMMINDEXER_HPP
#define DTK_COMMINDEXER_HPP

#include <map>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class CommIndexer
 * \brief Map the process ids of a local communicator into a global
 * communicator that encompasses it.
 */
//---------------------------------------------------------------------------//
class CommIndexer
{
  public:

    //@{
    //! Typedefs.
    typedef Teuchos::Comm<int>                             CommType;
    typedef Teuchos::RCP<const CommType>                   RCP_Comm;
    typedef std::map<int,int>                              IndexMap;
    //@}

  private:

    // Local to global process id map.
    IndexMap d_l2gmap;

  public:

    // Default constructor.
    CommIndexer();

    // Constructor.
    CommIndexer( RCP_Comm global_comm, RCP_Comm local_comm );

    // Destructor.
    ~CommIndexer();

    // Given a process id in the local communicator, return the distributed
    // object's process id in the global communicator.
    const int l2g( const int local_id ) const;

    //! Return the size of the local to global map.
    const int size() const
    { return d_l2gmap.size(); }
};

} // end namespace DataTransferKit

#endif // end DTK_COMMINDEXER_HPP

//---------------------------------------------------------------------------//
// end DTK_CommIndexer.hpp
//---------------------------------------------------------------------------//
