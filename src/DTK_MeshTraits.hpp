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
 * \class MeshTraits
 * \brief Mesh traits definitions.
 *
 * These traits correlate to the basic concept of a single topology mesh block
 * within DTK. A mesh block will consist of a globally unique list of vertex
 * ordinals of a type that implements Teuchos::OrdinalTraits (already
 * implemented for common ordinal types) and a set of globally unique element
 * ordinals of the same type. Vertices are described by a coordinate field
 * with coordinates of type double. Elements are described by a list of vertex
 * ordinals that designate their connectivity. For each element topology, the
 * order of the connecting vertices correlate to the canonical ordering scheme
 * for this particular mesh type. How this canonical ordering differs from DTK
 * canonical ordering is specified by the element connectivity permutation
 * list. Each block of mesh described by these traits must contain a single
 * element topology and a permutation list for that topology.
 */
//---------------------------------------------------------------------------//
template<typename MeshType>
class MeshTraits
{
  public:

    //@{
    //! Typedefs.
    //! Typedef for mesh type.
    typedef MeshType mesh_type;

    //! Typedef for global ordinal type. This type must implement
    //! Teuchos::OrdinalTraits.
    typedef typename MeshType::global_ordinal_type global_ordinal_type;

    //! Typedef for random access const iterator to vertex global_ordinal values.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const global_ordinal_type>
    const_vertex_iterator;

    //! Typedef for random access const iterator to coordinate values. This
    //! is enforcing a coordinate type of double.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const double>  
    const_coordinate_iterator;

    //! Typedef for random access const iterator to element global_ordinal
    //! values.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const global_ordinal_type>
    const_element_iterator;

    //! Typedef for random access const iterator to connectivity values.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const global_ordinal_type>
    const_connectivity_iterator;

    //! Typedef for random access const iterator to connectivity permutation
    //! list.
    typedef typename 
    std::iterator<std::random_access_iterator_tag, const int>
    const_permutation_iterator;
    //@}


    //@{
    //! Mesh vertex concepts.
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
    //! Mesh element concepts.
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
