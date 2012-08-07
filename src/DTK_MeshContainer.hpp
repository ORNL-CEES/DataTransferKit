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
 * \file DTK_MeshContainer.hpp
 * \author Stuart R. Slattery
 * \brief A simple default mesh container with a mesh traits
 * implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHCONTAINER_HPP
#define DTK_MESHCONTAINER_HPP

#include "DTK_MeshTraits.hpp"
#include "DTK_MeshTypes.hpp"

#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class MeshContainer
 * \brief A default mesh implementation.
 *
 * This container is used for rebuilding mesh data in the rendezvous
 * decomposition after serialization.
 */
//---------------------------------------------------------------------------//
template<typename GlobalOrdinal>
class MeshContainer
{

  public:

    //@{
    //! Typedefs.
    typedef GlobalOrdinal global_ordinal_type;
    //@}
    
    //! Default Constructor.
    MeshContainer()
    { /* ... */ }

    //! Constructor.
    MeshContainer( 
	const int vertex_dim,
	const Teuchos::ArrayRCP<GlobalOrdinal>& vertices,
	const Teuchos::ArrayRCP<const double>& coords,
	const DTK_ElementTopology element_topology,
	const int vertices_per_element,
	const Teuchos::ArrayRCP<GlobalOrdinal>& elements,
	const Teuchos::ArrayRCP<const GlobalOrdinal>& connectivity,
	const Teuchos::ArrayRCP<const int>& permutation_list )
	: d_vertex_dim( vertex_dim )
	, d_vertices( vertices )
	, d_coords( coords )
	, d_element_topology( element_topology )
	, d_vertices_per_element( vertices_per_element )
	, d_elements( elements )
	, d_connectivity( connectivity )
	, d_permutation_list( permutation_list )
    { /* ... */ }

    //! Destructor.
    ~MeshContainer()
    { /* ... */ }

    //! Get the dimension of the vertices.
    int getVertexDim() const
    { return d_vertex_dim; }

    //! Get the beginning of the vertices set.
    typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator 
    verticesBegin() const
    { return d_vertices.begin(); }

    //! Get the end of the vertices set.
    typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator 
    verticesEnd() const
    { return d_vertices.end(); }

    //! Get the beginning of the coordinates array.
    Teuchos::ArrayRCP<const double>::const_iterator coordsBegin() const
    { return d_coords.begin(); }

    //! Get the end of the coordinates array.
    Teuchos::ArrayRCP<const double>::const_iterator coordsEnd() const
    { return d_coords.end(); }

    //! Get the element topology.
    DTK_ElementTopology getElementTopology() const
    { return d_element_topology; }

    //! Get the number of vertices constructing a single element.
    int getVerticesPerElement() const
    { return d_vertices_per_element; }

    //! Get the beginning of the elements set.
    typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator 
    elementsBegin() const
    { return d_elements.begin(); }

    //! Get the end of the elements set.
    typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator 
    elementsEnd() const
    { return d_elements.end(); }

    //! Get the beginning of the connectivity array.
    typename Teuchos::ArrayRCP<const GlobalOrdinal>::const_iterator 
    connectivityBegin() const
    { return d_connectivity.begin(); }

    //! Get the ending of the connectivity array.
    typename Teuchos::ArrayRCP<const GlobalOrdinal>::const_iterator 
    connectivityEnd() const
    { return d_connectivity.end(); }

    //! Get the beginning of the permutation list.
    typename Teuchos::ArrayRCP<const int>::const_iterator 
    permutationBegin() const
    { return d_permutation_list.begin(); }

    //! Get the ending of the permutation list.
    typename Teuchos::ArrayRCP<const int>::const_iterator 
    permutationEnd() const
    { return d_permutation_list.end(); }
    
  private:

    // Vertex dimension.
    int d_vertex_dim;

    // Vertices.
    Teuchos::ArrayRCP<GlobalOrdinal> d_vertices;

    // Coordinates.
    Teuchos::ArrayRCP<const double> d_coords;

    // Element topology.
    DTK_ElementTopology d_element_topology;

    // Vertices per element.
    int d_vertices_per_element;

    // Elements.
    Teuchos::ArrayRCP<GlobalOrdinal> d_elements;

    // Connectivity.
    Teuchos::ArrayRCP<const GlobalOrdinal> d_connectivity;

    // Permutation list.
    Teuchos::ArrayRCP<const int> d_permutation_list;
};

//---------------------------------------------------------------------------//
/*
 * \class MeshTraits< MeshContainer<GlobalOrdinal> >
 * \brief MeshTraits specialization for the mesh container.
 */
//---------------------------------------------------------------------------//
template<>
template<typename GlobalOrdinal>
class MeshTraits< MeshContainer<GlobalOrdinal> >
{
  public:

    typedef MeshContainer<GlobalOrdinal> Container;
    typedef Container mesh_type;
    typedef typename Container::global_ordinal_type global_ordinal_type;

    typedef typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator 
    const_vertex_iterator;

    typedef Teuchos::ArrayRCP<const double>::const_iterator 
    const_coordinate_iterator;

    typedef typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator 
    const_element_iterator;

    typedef typename Teuchos::ArrayRCP<const GlobalOrdinal>::const_iterator 
    const_connectivity_iterator;

    typedef typename Teuchos::ArrayRCP<const int>::const_iterator 
    const_permutation_iterator;


    static inline int vertexDim( const Container& container )
    { return container.getVertexDim(); }

    static inline const_vertex_iterator 
    verticesBegin( const Container& container )
    { return container.verticesBegin(); }

    static inline const_vertex_iterator 
    verticesEnd( const Container& container )
    { return container.verticesEnd(); }

    static inline const_coordinate_iterator
    coordsBegin( const Container& container )
    { return container.coordsBegin(); }

    static inline const_coordinate_iterator
    coordsEnd( const Container& container )
    { return container.coordsEnd(); }


    static inline DTK_ElementTopology 
    elementTopology( const Container& container )
    { return container.getElementTopology(); }

    static inline int verticesPerElement( const Container& container )
    { return container.getVerticesPerElement(); }

    static inline const_element_iterator
    elementsBegin( const Container& container )
    { return container.elementsBegin(); }

    static inline const_element_iterator
    elementsEnd( const Container& container )
    { return container.elementsEnd(); }

    static inline const_connectivity_iterator
    connectivityBegin( const Container& container )
    { return container.connectivityBegin(); }

    static inline const_connectivity_iterator
    connectivityEnd( const Container& container )
    { return container.connectivityEnd(); }

    static inline const_permutation_iterator
    permutationBegin( const Container& container )
    { return container.permutationBegin(); }

    static inline const_permutation_iterator
    permutationEnd( const Container& container )
    { return container.permutationEnd(); }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MESHCONTAINER_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshContainer.hpp
//---------------------------------------------------------------------------//
