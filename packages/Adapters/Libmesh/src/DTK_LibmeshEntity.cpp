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
 * \brief DTK_LibmeshEntity.cpp
 * \author Stuart R. Slattery
 * \brief Libmesh entity interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_LibmeshEntity.hpp"
#include "DTK_LibmeshEntityImpl.hpp"

#include <libmesh/elem.h>
#include <libmesh/node.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor. Elem specialization.
template <>
LibmeshEntity<libMesh::Elem>::LibmeshEntity(
    const Teuchos::Ptr<libMesh::Elem> &libmesh_elem,
    const Teuchos::Ptr<libMesh::MeshBase> &libmesh_mesh,
    const Teuchos::Ptr<LibmeshAdjacencies> &adjacencies )
{
    this->b_entity_impl = Teuchos::rcp( new LibmeshEntityImpl<libMesh::Elem>(
        libmesh_elem, libmesh_mesh, adjacencies ) );
}

//---------------------------------------------------------------------------//
// Constructor. Node specialization.
template <>
LibmeshEntity<libMesh::Node>::LibmeshEntity(
    const Teuchos::Ptr<libMesh::Node> &libmesh_node,
    const Teuchos::Ptr<libMesh::MeshBase> &libmesh_mesh,
    const Teuchos::Ptr<LibmeshAdjacencies> &adjacencies )
{
    this->b_entity_impl = Teuchos::rcp( new LibmeshEntityImpl<libMesh::Node>(
        libmesh_node, libmesh_mesh, adjacencies ) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_Entity.cpp
//---------------------------------------------------------------------------//
