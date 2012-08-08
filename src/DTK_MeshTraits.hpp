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
 * \brief DTK_MeshTraits.hpp
 * \author Stuart R. Slattery
 * \brief Declaration of mesh traits.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHTRAITS_HPP
#define DTK_MESHTRAITS_HPP

#include <iterator>

#include "DTK_MeshTypes.hpp"

namespace DataTransferKit
{

/*!
 * \brief Dummy struct. If a type does not create a specialization this will
 * not compile.
 */
template<typename UndefinedMeshType>
struct UndefinedMeshTraits
{
    static inline UndefinedMeshType notDefined() 
    { return UndefinedMeshType::this_type_is_missing_a_specialization(); }
};

//---------------------------------------------------------------------------//
/*!
  \class MeshTraits
  \brief Mesh traits definitions.
 
  In order to access DTK mesh services, a subset of the information needed to
  describe the mesh is required. This subset consists of vertices and their
  coordinates, elements and the vertices that construct them, and the
  communicator over which they are defined. The vertices that construct an
  element have both a canonical ordering consistent across all elements of
  that topology in a mesh and a permutation list that describes how this
  ordering varies from DTK canonical ordering. 

  Vertices are the lowest level geometric component of the mesh. All vertices
  have a globally unique ordinal serving as an identification number for the
  vertex in global operations. A vertex can have 1, 2, or 3 dimensions but all
  vertices in a mesh must have the same dimension. To specify its geometric
  position, each vertex has Cartesian (x,y,z) coordinates. A vertex must
  provide only the coordinates for the specified vertex dimension, no more or
  no less (e.g. a 2 dimensional vertex must provide x and y coordinates but
  not a z coordinate). A vertex may be repeated any number of times across the
  parallel domain with unlimited local and global instances. However, every
  vertex with the same globally unique ordinal must have the same
  coordinates. We make a distinction here between vertices and nodes. In the
  context of DTK, a vertex is purely a geometric object. It describes the
  spatial positioning and geometric bounds of an element. A node is purely a
  mathematical object. It descrbibes the descretization associated with a
  particular element described within the natural coordinate system of that
  element. It is possible that in the physical coordinate frame that a node
  and vertex may occupy the same geometric location, however DTK does not
  consider nodes in its formulation.

  Elements are the second level of abstraction in the mesh description above
  vertices. All elements have a globally unique ordinal serving as an
  identification number for the element in global operations. This globally
  unique ordinal can be the same as a globally unique ordinal for a vertex in
  the mesh as DTK distinguishes between vertices and elements. An element has
  a topology defining its physical structure (e.g. tetrahedron, hexahedron,
  etc.) and a number of vertices needed to generate that topology. Elements
  are constructed from vertices via a connectivity list. The connectivity list
  for a particular element will contain the unique vertex global ordinals that
  construct its linear form. An element may be repeated any number of times
  across the parallel domain (i.e. it may have unlimited local instances),
  however, every globally unique ordinal must have the same connectivity list
  associated with it. For consistency, DTK uses the MoaB Canonical Numbering
  (MBCN) scheme as a canonical ordering scheme. Each element in a client mesh
  can be described with a connectivity list using any canonical scheme of
  choice, however, the relationship between this canonical numbering scheme
  and the DTK canonical numbering scheme must be made available. Each element
  topology is therefore also described by a permutation list. A permutation
  list specifies the variation in ordering between the DTK canonical numbering
  scheme and the client canonical numbering scheme. A permutation list must be
  described globally, regardless of whether or not elements exist on a
  particular process. See DTK_ElementTopology for canonical element topologies
  as defined by DTK. Mesh elements may not intersect any other elements in a
  single mesh description. An element may intersect other elements if those
  elements exist in another mesh (this is in fact a common situation in data
  transfer).

  MeshTraits correlate to the basic concept of a single topology mesh block
  within DTK. They have the following properties:

  Mesh vertices have D dimensions and may not exceed three dimensions \f$
  \Big\{ d_0, ..., d_D \Big\} \f$. Vertices are identified by a unique global
  ordinal. If there are N vertices in the mesh then their ordinals are given
  as \f$ \Big\{ n^0, n^1, n^2, ..., n^N \Big\} \f$. Vertex coordinates are
  blocked by dimension such that if there are N vertices in the mesh block
  then they are stored as: 

  \f[ 
  \Big\{ x^0_0, x^1_0, x^2_0, ..., x^N_0, x^0_1, x^1_1, x^2_1, ... x^N_1,
  x^0_2, x^1_2, x^2_2, ... x^N_2 \Big\}
  \f]

  with the superscript denoting which vertex the coordinate corresponds to and
  the subscript denoting which dimension the coordinate is for.  The ordering
  of the vertices is implicilty bound to the global ordinals such that the
  coordinates for a vertex with ordinal \f$ n_N \f$ are \f$ \Big\{ x^N_0,
  x^N_1, x^N_2 \Big\} \f$.

  Mesh elements have a topology defined by a DTK_ElementTopology enumeration
  with a specified number of vertices, P, needed to construct the
  topology. Elements are identified by a unique global ordinal. If there are M
  elements in the mesh then their ordinals are given as \f$ \Big\{ m^0, m^1,
  m^2, ..., m^M \Big\} \f$. The connecting vertices for the elements in the
  mesh are defined using the vertex global ordinals such that an element, m,
  can be described with a list \f$ \Big\{ n_0^m, n_1^m, ..., n_P^m \Big\}
  \f$. The connectivity information is accessed by blocks in the same manner
  as coordinates such that: 

  \f[ 
  \Big\{ n_0^0, n_0^1, n_0^2, ..., n_0^M, n_1^0, n_1^1, n_1^2, ..., n_1^M,
  ..., n_P^0, n_P^1, ..., n_P^M \Big\} 
  \f]
 
  describes the connectivity of a mesh block with the superscript
  cooresponding to which element the vertex ordinal constructs and the
  superscript corresponding to the canonical vertex index for the given
  element topology.  Finally, a permutation list defines the difference in
  ordering between a client element topology connectivity ordering and DTK
  canonical ordering. This list, defined as \f$ \Big\{ p_0, p_1, ..., p_P
  \Big\}\f$, must be defined for every instance of the mesh, regardless of
  whether or not the mesh contains any data. Here, the entry \f$ p_P \f$ gives
  which canonical vertex index in the client connectivity list cooresponds to
  the \f$ P^{th} \f$ vertex in the DTK canonical vertex list for that
  topology.
*/
//---------------------------------------------------------------------------//
template<typename MeshType>
class MeshTraits
{
  public:

    //@{
    //! Typedef for mesh type.
    typedef MeshType mesh_type;

    //! Typedef for global ordinal type. This type must implement
    //  Teuchos::OrdinalTraits.
    typedef typename MeshType::global_ordinal_type global_ordinal_type;

    //! Typedef for random access const iterator to vertex global ordinal
    //  values.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const global_ordinal_type>
    const_vertex_iterator;

    //! Typedef for random access const iterator to coordinate
    //  values. Coordinates are required to be of type double.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const double>  
    const_coordinate_iterator;

    //! Typedef for random access const iterator to element global ordinal
    //  values.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const global_ordinal_type>
    const_element_iterator;

    //! Typedef for random access const iterator to connectivity values.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const global_ordinal_type>
    const_connectivity_iterator;

    //! Typedef for random access const iterator to connectivity permutation
    // list.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const int>
    const_permutation_iterator;
    //@}


    //@{
    // Mesh vertex concepts.
    /*!
     * \brief Return the dimension of the vertices in this mesh block.
     */
    static inline int vertexDim( const MeshType& mesh_block )
    { UndefinedMeshTraits<MeshType>::notDefined(); return 0; }

    /*!
     * \brief Return the const iterator to the beginning of the vertex global
     * ordinal block in this mesh block.
     */
    static inline const_vertex_iterator 
    verticesBegin( const MeshType& mesh_block )
    { UndefinedMeshTraits<MeshType>::notDefined(); return 0; }

    /*!
     * \brief Return the const iterator to the end of the vertex global ordinal
     * block in this mesh block.
     */ 
    static inline const_vertex_iterator 
    verticesEnd( const MeshType& mesh_block )
    { UndefinedMeshTraits<MeshType>::notDefined(); return 0; }

    /*!
     * \brief Return the const iterator to the beginning of the vertex
     * coordinate block in this mesh block. These coordinates are required to
     * be three dimensional and blocked.
     * ( x0, x1, x2, ... , xN, y0, y1, y2, ... , yN, z0, z1, z2, ... , zN )
     */
    static inline const_coordinate_iterator 
    coordsBegin( const MeshType& mesh_block )
    { UndefinedMeshTraits<MeshType>::notDefined(); return 0; }

    /*!
     * \brief Return the const iterator to the end of the vertex coordinate
     * block in this mesh block. These coordinates are requried to be three
     * dimensional and blocked.
     * ( x0, x1, x2, ... , xN, y0, y1, y2, ... , yN, z0, z1, z2, ... , zN )
     */
    static inline const_coordinate_iterator 
    coordsEnd( const MeshType& mesh_block )
    { UndefinedMeshTraits<MeshType>::notDefined(); return 0; }
    //@}


    //@{
    // Mesh element concepts.
    /*! 
     * \brief Return the element topology for this mesh block
     * (DTK_ElementTopology enum).
     */
    static inline DTK_ElementTopology 
    elementTopology( const MeshType& mesh_block )
    { UndefinedMeshTraits<MeshType>::notDefined(); return 0;}

    /*! 
     * \brief Return the number of vertices that constructs an individual
     * element in this mesh block. All elements in the mesh must be
     * constructed with the same number of vertices.
     */
    static inline int verticesPerElement( const MeshType& mesh_block )
    { UndefinedMeshTraits<MeshType>::notDefined(); return 0; }

    /*! 
     * \brief Return the const iterator to the beginning of the element global
     * ordinal block in this mesh block.
     */
    static inline const_element_iterator 
    elementsBegin( const MeshType& mesh_block )
    { UndefinedMeshTraits<MeshType>::notDefined(); return 0; }

    /*! 
     * \brief Return the const iterator to the end of the element global
     * ordinal block in this mesh block.
     */
    static inline const_element_iterator 
    elementsEnd( const MeshType& mesh_block )
    { UndefinedMeshTraits<MeshType>::notDefined(); return 0; }

    /*! 
     * \brief Return the const iterator to the beginning of the element
     * connectivity block in this mesh block. The connectivity entries are
     * required to be blocked. 
     * ( element0( c0 ), element1( c0 ), ... , elementN( c0 ), element0( c1 ),
     * element1( c1 ), ... , elementN( c1 ), ... , elementN( cn ) )
     */
    static inline const_connectivity_iterator 
    connectivityBegin( const MeshType& mesh_block )
    { UndefinedMeshTraits<MeshType>::notDefined(); return 0; }

    /*! 
     * \brief Return the const iterator to the end of the element connectivity
     * block in this mesh block. The connectivity entries are required to be
     * blocked.  
     * ( element0( c0 ), element1( c0 ), ... , elementN( c0 ), element0( c1 ),
     * element1( c1 ), ... , elementN( c1 ), ... , elementN( cn ) )
     */
    static inline const_connectivity_iterator 
    connectivityEnd( const MeshType& mesh_block )
    { UndefinedMeshTraits<MeshType>::notDefined(); return 0; }

    /*! 
     * \brief Return the const iterator to the beginning of the element
     * connectivity permutation list.
     */
    static inline const_permutation_iterator
    permutationBegin( const MeshType& mesh_block )
    { UndefinedMeshTraits<MeshType>::notDefined(); return 0; }

    /*! 
     * \brief Return the const iterator to the end of the element connectivity
     * permutation list. 
     */
    static inline const_permutation_iterator
    permutationEnd( const MeshType& mesh_block )
    { UndefinedMeshTraits<MeshType>::notDefined(); return 0; }
    //@}
};

} // end namespace DataTransferKit

#endif // end DTK_MESHTRAITS_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshTraits.hpp
//---------------------------------------------------------------------------//
