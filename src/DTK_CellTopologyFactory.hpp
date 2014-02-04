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
 * \file DTK_CellTopologyFactory.hpp
 * \author Stuart R. Slattery
 * \brief Factory method declaration for cell topology data.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CELLTOPOLOGYFACTORY_HPP
#define DTK_CELLTOPOLOGYFACTORY_HPP

#include <DTK_MeshTypes.hpp>

#include <Teuchos_RCP.hpp>

#include <Shards_CellTopology.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class CellTopologyFactory 
 * \brief A stateless factory class for creating Shards::CellTopology
 * structures based on a moab::EntityType.
 */
//---------------------------------------------------------------------------//
class CellTopologyFactory
{
  public:

    //@{
    //! Typedefs.
    typedef Teuchos::RCP<shards::CellTopology>   RCP_CellTopology;
    //@}

    // Constructor.
    CellTopologyFactory();

    // Destructor.
    ~CellTopologyFactory();

    // Factory method.
    static RCP_CellTopology create( const DTK_ElementTopology element_topology, 
				    const int num_element_vertices );
};

} // end namespace DataTransferKit

#endif // end DTK_CELLTOPOLOGYFACTORY_HPP

//---------------------------------------------------------------------------//
// end DTK_CellTopologyFactory.hpp
//---------------------------------------------------------------------------//
