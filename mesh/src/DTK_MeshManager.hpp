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
 * \file DTK_MeshManager.hpp
 * \author Stuart R. Slattery
 * \brief Mesh manager declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHMANAGER_HPP
#define DTK_MESHMANAGER_HPP

#include "DTK_MeshTraits.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{

template<class Mesh>
class MeshManager
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh                                                mesh_type;
    typedef MeshTraits<Mesh>                                    MT;
    typedef Teuchos::Comm<int>                                  CommType;
    typedef Teuchos::RCP<const CommType>                        RCP_Comm;
    typedef typename Teuchos::ArrayRCP<Mesh>::const_iterator    BlockIterator;
    //@}

    // Constructor.
    MeshManager( const Teuchos::ArrayRCP<Mesh>& mesh_blocks,
		 const RCP_Comm& comm,
		 const std::size_t dim );

    // Destructor.
    ~MeshManager();

    // Get the number of mesh blocks.
    int getNumBlocks() const
    { return d_mesh_blocks.size(); }

    // Get the iterator to the beginning of the mesh blocks.
    BlockIterator blocksBegin() const
    { return d_mesh_blocks.begin(); }

    // Get the iterator to the end of the mesh blocks.
    BlockIterator blocksEnd() const
    { return d_mesh_blocks.end(); }

    // Get the communicator for the mesh.
    const RCP_Comm& getComm() const
    { return d_comm; }

    // Get the physical dimension of the mesh.
    std::size_t getDim() const
    { return d_dim; }

  private:

    // Validate the mesh.
    void validate();

  private:

    // Mesh block array.
    Teuchos::ArrayRCP<Mesh> d_mesh_blocks;

    // Communicator over which the mesh is defined.
    RCP_Comm d_comm;

    // The physical dimension of the mesh.
    std::size_t d_dim;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_MeshManager_def.hpp"

#endif // end DTK_MESHMANAGER_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshManager.hpp
//---------------------------------------------------------------------------//

