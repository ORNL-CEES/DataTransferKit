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

#include <boost/tr1/unordered_map.hpp>

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
 \class MeshManager
 \brief Manager object for mesh.
 
 The MeshManager manages the topology blocks that construct a mesh and its
 parallel decomposition. The mesh manager also keeps track of the mesh
 vertices and elements that are active in partitioning and searching. Those
 that are not active are typically outside some domain of interest.

 All elements in a block must have the same topology and number of vertices. A
 mesh may contain as many blocks as desired. Multiple blocks with the same
 mesh topology may exist. Vertices may be repeated in different mesh blocks
 provided that they maintain the same unique global ordinal and
 coordinates. Elements may be repeated in different mesh blocks provided that
 they maintain the same unique global ordinal and connectivity list. All
 elements and vertices in all blocks of the mesh must have the same
 dimension. Multiple mesh blocks may exist in the same spatial region as they
 are merely a means of subdividing the mesh into groups of elements based on
 their topology. This behavior will be the common when hybrid meshes are
 involved (e.g. a mesh that contains hexahedrons and tetrahedrons). Mesh
 blocks may be either structured or unstructured.

 Mesh blocks and the elements and vertices they contain may be partitioned in
 any fashion provided that all vertices, elements, and blocks of a mesh
 description exist in a communication space operated on by the same parallel
 communicator. Different blocks in a single mesh description may not have
 different communicators. Each mesh description may have its own
 communicator. Global knowledge of the parallel decomposition of a given mesh
 description is not required. Only local mesh data access along with the
 proper communicator is required.
 */
//---------------------------------------------------------------------------//
template<class Mesh>
class MeshManager
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh                                                mesh_type;
    typedef Teuchos::RCP<Mesh>                                  RCP_Mesh;
    typedef MeshTraits<Mesh>                                    MT;
    typedef typename MT::global_ordinal_type                    GlobalOrdinal;
    typedef Teuchos::Comm<int>                                  CommType;
    typedef Teuchos::RCP<const CommType>                        RCP_Comm;
    typedef Teuchos::ArrayRCP<RCP_Mesh>                         BlockArray;
    typedef typename BlockArray::const_iterator                 BlockIterator;
    //@}

    // Constructor.
    MeshManager( const Teuchos::ArrayRCP<RCP_Mesh>& mesh_blocks,
		 const RCP_Comm& comm,
		 const int dim );

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

    //! Get a block of mesh.
    const RCP_Mesh& getBlock( const int block_id ) const
    { return d_mesh_blocks[ block_id ]; }

    //! Get the communicator for the mesh.
    const RCP_Comm& comm() const
    { return d_comm; }

    //! Get the physical dimension of the mesh.
    int dim() const
    { return d_dim; }

    //! Set the active vertices for a block.
    void setActiveVertices( const Teuchos::Array<short int>& active_vertices,
			    const int block_id )
    { d_active_vertices[ block_id ] = active_vertices; }

    //! Set the active elements for a block.
    void setActiveElements( const Teuchos::Array<short int>& active_elements,
			    const int block_id )
    { d_active_elements[ block_id ] = active_elements; }

    //! Get the active vertices for a block.
    Teuchos::ArrayView<short int> getActiveVertices( const int block_id )
    { return d_active_vertices[ block_id ](); }

    //! Get the active elements for a block.
    Teuchos::ArrayView<short int> getActiveElements( const int block_id )
    { return d_active_elements[ block_id ](); }

    // Compute the global bounding box around the entire mesh.
    BoundingBox globalBoundingBox();

    // Building indexing for faster coordinate access via connectivity.
    void buildIndexing();

    // Given the local id of an element in the mesh, get its block id and its
    // local id in the block.
    void getLocalElementIds( const int local_elem_id,
			     int& block_id,
			     int& block_elem_id ) const;

    // Given the block id of a vertex and its global id, get its local id in
    // that block.
    int getLocalVertexId( const int block_id, const int vertex_gid ) const;

  private:

    // Validate the mesh to the domain model.
    void validate();

  private:

    // Mesh block array.
    Teuchos::ArrayRCP<RCP_Mesh> d_mesh_blocks;

    // Communicator over which the mesh is defined.
    RCP_Comm d_comm;

    // The physical dimension of the mesh.
    int d_dim;

    // Active vertices in each mesh block.
    Teuchos::Array< Teuchos::Array<short int> > d_active_vertices;

    // Active elements in each mesh block.
    Teuchos::Array< Teuchos::Array<short int> > d_active_elements;

    // Array of cumulative number elements in each block.
    Teuchos::Array<GlobalOrdinal> d_cumulative_elements;

    // Vertex global-to-local indexers for each block.
    Teuchos::Array<std::tr1::unordered_map<GlobalOrdinal,int> > d_vertex_g2l;
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

