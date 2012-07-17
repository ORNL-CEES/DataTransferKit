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
#include "DTK_BoundingBox.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class MeshManager
 * \brief Manager object for mesh.
 *
 * The mesh manager manages the blocks that construct a mesh and its parallel
 * decomposition.
 */
//---------------------------------------------------------------------------//
template<class Mesh>
class MeshManager
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh                                                mesh_type;
    typedef MeshTraits<Mesh>                                    MT;
    typedef typename MT::global_ordinal_type                    GlobalOrdinal;
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

    //! Get the number of mesh blocks.
    int getNumBlocks() const
    { return d_mesh_blocks.size(); }

    // Get the local number of elements in the mesh.    
    GlobalOrdinal localNumElements() const;

    // Get the global number of elements in the mesh.    
    GlobalOrdinal globalNumElements() const;

    //! Get the iterator to the beginning of the mesh blocks.
    BlockIterator blocksBegin() const
    { return d_mesh_blocks.begin(); }

    //! Get the iterator to the end of the mesh blocks.
    BlockIterator blocksEnd() const
    { return d_mesh_blocks.end(); }

    //! Get the communicator for the mesh.
    const RCP_Comm& comm() const
    { return d_comm; }

    //! Get the physical dimension of the mesh.
    std::size_t dim() const
    { return d_dim; }

    //! Set the active nodes for a block.
    void setActiveNodes( const Teuchos::Array<short int>& active_nodes,
			 const int block_id )
    { d_active_nodes[ block_id ] = active_nodes; }

    //! Set the active elements for a block.
    void setActiveElements( const Teuchos::Array<short int>& active_elements,
			    const int block_id )
    { d_active_elements[ block_id ] = active_elements; }

    //! Get the active nodes for a block.
    Teuchos::ArrayView<short int> getActiveNodes( const int block_id )
    { return d_active_nodes[ block_id ](); }

    //! Get the active elements for a block.
    Teuchos::ArrayView<short int> getActiveElements( const int block_id )
    { return d_active_elements[ block_id ](); }

    // Compute the global bounding box around the entire mesh.
    BoundingBox globalBoundingBox();

  private:

    // Validate the mesh to the domain model.
    void validate();

  private:

    // Mesh block array.
    Teuchos::ArrayRCP<Mesh> d_mesh_blocks;

    // Communicator over which the mesh is defined.
    RCP_Comm d_comm;

    // The physical dimension of the mesh.
    std::size_t d_dim;

    // Active nodes in each mesh block.
    Teuchos::Array< Teuchos::Array<short int> > d_active_nodes;

    // Active elements in each mesh block.
    Teuchos::Array< Teuchos::Array<short int> > d_active_elements;
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

