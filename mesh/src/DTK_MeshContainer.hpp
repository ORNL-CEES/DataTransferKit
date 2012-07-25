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
 * \brief A simple mesh container for rebuilding mesh data after
 * serialization.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHCONTAINER_HPP
#define DTK_MESHCONTAINER_HPP

#include <iterator>

#include "DTK_MeshTraits.hpp"

#include <Teuchos_ArrayRCP.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class MeshContainer
 * \brief A container for rebuilding mesh data after serialization.
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
    MeshContainer( const int node_dim,
		   const Teuchos::ArrayRCP<GlobalOrdinal>& nodes,
		   const Teuchos::ArrayRCP<const double>& coords,
		   const int element_topology,
		   const int nodes_per_element,
		   const Teuchos::ArrayRCP<GlobalOrdinal>& elements,
		   const Teuchos::ArrayRCP<const GlobalOrdinal>& connectivity,
		   const Teuchos::ArrayRCP<const std::size_t>& permutation_list )
	: d_node_dim( node_dim )
	, d_nodes( nodes )
	, d_coords( coords )
	, d_element_topology( element_topology )
	, d_nodes_per_element( nodes_per_element )
	, d_elements( elements )
	, d_connectivity( connectivity )
	, d_permutation_list( permutation_list )
    { /* ... */ }

    //! Destructor.
    ~MeshContainer()
    { /* ... */ }

    //! Get the dimension of the nodes.
    std::size_t getNodeDim() const
    { return d_node_dim; }

    //! Get the beginning of the nodes set.
    typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator 
    nodesBegin() const
    { return d_nodes.begin(); }

    //! Get the end of the nodes set.
    typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator 
    nodesEnd() const
    { return d_nodes.end(); }

    //! Get the beginning of the coordinates array.
    Teuchos::ArrayRCP<const double>::const_iterator coordsBegin() const
    { return d_coords.begin(); }

    //! Get the end of the coordinates array.
    Teuchos::ArrayRCP<const double>::const_iterator coordsEnd() const
    { return d_coords.end(); }

    //! Get the element topology.
    std::size_t getElementTopology() const
    { return d_element_topology; }

    //! Get the number of nodes constructing a single element.
    std::size_t getNodesPerElement() const
    { return d_nodes_per_element; }

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
    typename Teuchos::ArrayRCP<const std::size_t>::const_iterator 
    permutationBegin() const
    { return d_permutation_list.begin(); }

    //! Get the ending of the permutation list.
    typename Teuchos::ArrayRCP<const std::size_t>::const_iterator 
    permutationEnd() const
    { return d_permutation_list.end(); }
    
  private:

    // Node dimension.
    std::size_t d_node_dim;

    // Nodes.
    Teuchos::ArrayRCP<GlobalOrdinal> d_nodes;

    // Coordinates.
    Teuchos::ArrayRCP<const double> d_coords;

    // Element topology.
    std::size_t d_element_topology;

    // Nodes per element.
    std::size_t d_nodes_per_element;

    // Elements.
    Teuchos::ArrayRCP<GlobalOrdinal> d_elements;

    // Connectivity.
    Teuchos::ArrayRCP<const GlobalOrdinal> d_connectivity;

    // Permutation list.
    Teuchos::ArrayRCP<const std::size_t> d_permutation_list;
};

//---------------------------------------------------------------------------//
/*
 * \class MeshTraits
 * \brief MeshTraits specialization for the mesh container.
 */
//---------------------------------------------------------------------------//
template<>
template<typename GlobalOrdinal>
class MeshTraits< MeshContainer<GlobalOrdinal> >
{
    typedef MeshContainer<GlobalOrdinal> Container;
    typedef typename Container::global_ordinal_type global_ordinal_type;

    typedef typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator 
    const_node_iterator;

    typedef Teuchos::ArrayRCP<const double>::const_iterator 
    const_coordinate_iterator;

    typedef typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator 
    const_element_iterator;

    typedef typename Teuchos::ArrayRCP<const GlobalOrdinal>::const_iterator 
    const_connectivity_iterator;

    typedef typename Teuchos::ArrayRCP<const std::size_t>::const_iterator 
    const_permutation_iterator;


    static inline std::size_t nodeDim( const Container& container )
    { return container.getNodeDim(); }

    static inline const_node_iterator nodesBegin( const Container& container )
    { return container.nodesBegin(); }

    static inline const_node_iterator 
    nodesEnd( const Container& container )
    { return container.nodesEnd(); }

    static inline const_coordinate_iterator
    coordsBegin( const Container& container )
    { return container.coordsBegin(); }

    static inline const_coordinate_iterator
    coordsEnd( const Container& container )
    { return container.coordsEnd(); }


    static inline std::size_t elementTopology( const Container& container )
    { return container.getElementTopology(); }

    static inline std::size_t nodesPerElement( const Container& container )
    { return container.getNodesPerElement(); }

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
