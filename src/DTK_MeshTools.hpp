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
 * \file DTK_MeshTools.hpp
 * \author Stuart R. Slattery
 * \brief MeshTools declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHTOOLS_HPP
#define DTK_MESHTOOLS_HPP

#include "DTK_MeshTypes.hpp"
#include "DTK_MeshBlock.hpp"
#include "DTK_BoundingBox.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class MeshTools
 * \brief A stateless class with tools for operating on single mesh blocks.
 */ 
//---------------------------------------------------------------------------//
class MeshTools
{
  public:

    //@{
    //! Typedefs.
    typedef Teuchos::Comm<int>                  CommType;
    typedef Teuchos::RCP<const CommType>        RCP_Comm;
    //@}

    //! Constructor.
    MeshTools()
    { /* ... */ }

    //! Destructor.
    ~MeshTools()
    { /* ... */ }

    // Get the local bounding box for a mesh block.
    static BoundingBox localBoundingBox( const Teuchos::RCP<MeshBlock>& mesh );

    // Get the global bounding box for a mesh block over the given
    // communicator. 
    static BoundingBox globalBoundingBox( const Teuchos::RCP<MeshBlock>& mesh, 
					  const RCP_Comm& comm );
};

//---------------------------------------------------------------------------//

} // end namepsace DataTransferKit

#endif // end DTK_MESHTOOLS_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshTools.hpp
//---------------------------------------------------------------------------//
