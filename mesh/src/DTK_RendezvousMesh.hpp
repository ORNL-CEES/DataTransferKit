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
 * \file DTK_RendezvousMesh.hpp
 * \author Stuart R. Slattery
 * \brief Concrete rendezvous mesh declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_RENDEZVOUSMESH_HPP
#define DTK_RENDEZVOUSMESH_HPP

#include <map>

#include "DTK_MeshTraits.hpp"
#include "DTK_MeshManager.hpp"

#include <MBInterface.hpp>
#include <MBRange.hpp>

#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class RendezvousMesh
 * \brief Concrete mesh data structure for building mesh in the rendezvous
 * decomposition.
 */
//---------------------------------------------------------------------------//
template<typename GlobalOrdinal>
class RendezvousMesh
{
  public:

    //@{
    //! Typedefs.
    typedef GlobalOrdinal                               global_ordinal_type;
    typedef Teuchos::RCP<moab::Interface>               RCP_Moab;
    typedef std::map<moab::EntityHandle,GlobalOrdinal>  HandleMap;
    //@}

    // Constructor.
    RendezvousMesh( const RCP_Moab& moab, 
		    const moab::Range& elements,
		    const HandleMap& handle_map );

    // Destructor.
    ~RendezvousMesh();

    //! Get the Moab interface.
    const RCP_Moab& getMoab() const
    { return d_moab; }

    //! Get the mesh elements.
    const moab::Range& getElements() const
    { return d_elements; }

    //! Given a moab element handle return the corresponding native element
    //! global ordinal.
    GlobalOrdinal getNativeOrdinal( const moab::EntityHandle& moab_handle ) const
    { return d_handle_map.find( moab_handle )->second; }

  private:

    //! Moab interface implementation.
    RCP_Moab d_moab;

    //! Moab mesh elements.
    moab::Range d_elements;

    //! Moab element handle to native element handle map.
    HandleMap d_handle_map;
};

//---------------------------------------------------------------------------//
//! DTK to MOAB topology translation table.
const moab::EntityType moab_topology_table[] =
{
    moab::MBVERTEX,  // DTK_VERTEX
    moab::MBEDGE,    // DTK_LINE_SEGMENT
    moab::MBTRI,     // DTK_TRIANGLE
    moab::MBQUAD,    // DTK_QUADRILATERAL
    moab::MBTET,     // DTK_TETRAHEDRON
    moab::MBPYRAMID, // DTK_PYRAMID
    moab::MBPRISM,   // DTK_WEDGE
    moab::MBHEX      // DTK_HEXAHEDRON
};

//---------------------------------------------------------------------------//
//! MOAB ErrorCode table.
const std::string moab_error_table[] =
{
    "MB_SUCCESS",
    "MB_INDEX_OUT_OF_RANGE",
    "MB_TYPE_OUT_OF_RANGE",
    "MB_MEMORY_ALLOCATION_FAILED",
    "MB_ENTITY_NOT_FOUND",
    "MB_MULTIPLE_ENTITIES_FOUND",
    "MB_TAG_NOT_FOUND",
    "MB_FILE_DOES_NOT_EXIST",
    "MB_FILE_WRITE_ERROR",
    "MB_NOT_IMPLEMENTED",
    "MB_ALREADY_ALLOCATED",
    "MB_VARIABLE_DATA_LENGTH",
    "MB_INVALID_SIZE",
    "MB_UNSUPPORTED_OPERATION",
    "MB_UNHANDLED_OPTION",
    "MB_STRUCTURED_MESH",
    "MB_FAILURE"
};

//---------------------------------------------------------------------------//
// Non-member creation methods.
//---------------------------------------------------------------------------//

// Create a RendezvousMesh from mesh manager.
template<class Mesh>
Teuchos::RCP< RendezvousMesh<typename MeshTraits<Mesh>::global_ordinal_type> >
createRendezvousMesh( const MeshManager<Mesh>& mesh_manager );

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_RendezvousMesh_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_RENDEZVOUSMESH_HPP

//---------------------------------------------------------------------------//
// end DTK_RendezvousMesh.hpp
//---------------------------------------------------------------------------//

