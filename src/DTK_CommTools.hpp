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
 * \file DTK_CommTools.hpp
 * \author Stuart R. Slattery
 * \brief CommTools declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_COMMTOOLS_HPP
#define DTK_COMMTOOLS_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class CommTools
 * \brief A stateless class with tools for operating on communicators.
 */ 
//---------------------------------------------------------------------------//
class CommTools
{
  public:

    //@{
    //! Typedefs.
    typedef Teuchos::Comm<int>                  CommType;
    typedef Teuchos::RCP<const CommType>        RCP_Comm;
    //@}

    //! Constructor.
    CommTools()
    { /* ... */ }

    //! Destructor.
    ~CommTools()
    { /* ... */ }

    // Get comm world.
    static void getCommWorld( RCP_Comm& comm_world );

    // Check whether two communicators own the same communication space.
    static bool equal( const RCP_Comm& comm_A, 
                       const RCP_Comm& comm_B,
                       const RCP_Comm& comm_global = Teuchos::null );

    // Generate the union of two communicators.
    static void unite( const RCP_Comm& comm_A, 
                       const RCP_Comm& comm_B,
		       RCP_Comm& comm_union,
                       const RCP_Comm& comm_global = Teuchos::null );

    // Generate the intersection of two communicators.
    static void intersect( const RCP_Comm& comm_A, 
                           const RCP_Comm& comm_B,
			   RCP_Comm& comm_intersection,
                           const RCP_Comm& comm_global = Teuchos::null );
};

//---------------------------------------------------------------------------//

} // end namepsace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_COMMTOOLS_HPP

//---------------------------------------------------------------------------//
// end DTK_CommTools.hpp
//---------------------------------------------------------------------------//
