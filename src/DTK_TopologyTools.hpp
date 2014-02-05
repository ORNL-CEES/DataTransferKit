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
 * \file DTK_TopologyTools.hpp
 * \author Stuart R. Slattery
 * \brief TopologyTools declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_TOPOLOGYTOOLS_HPP
#define DTK_TOPOLOGYTOOLS_HPP

#include <DTK_BoundingBox.hpp>

#include <MBInterface.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class TopologyTools
 * \brief A stateless class with tools for operating on mesh topologies.
 */ 
//---------------------------------------------------------------------------//
class TopologyTools
{
  public:

    //! Constructor.
    TopologyTools()
    { /* ... */ }

    //! Destructor.
    ~TopologyTools()
    { /* ... */ }

    // Get the center of the reference cell of the given topology.
    template<typename MDArray>
    static void referenceCellCenter( const shards::CellTopology& cell_topo,
				     MDArray& center );

    // Box-element overlap query.
    static bool boxElementOverlap( const BoundingBox& box,
				   const moab::EntityHandle element,
				   const Teuchos::RCP<moab::Interface>& moab );

    // Element-in-geometry query.
    template<class Geometry>
    static bool elementInGeometry( const Geometry& geometry,
				   const moab::EntityHandle element,
				   const Teuchos::RCP<moab::Interface>& moab,
				   const double tolerance,
				   bool all_vertices_for_inclusion );
};

} // end namepsace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_TopologyTools_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_TOPOLOGYTOOLS_HPP

//---------------------------------------------------------------------------//
// end DTK_TopologyTools.hpp
//---------------------------------------------------------------------------//
