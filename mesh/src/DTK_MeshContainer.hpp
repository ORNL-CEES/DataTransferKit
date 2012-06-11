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

template<typename Handle>
class MeshContainer
{

  public:

    //@{
    //! Typedefs.
    typedef Handle handle_type;
    //@}
    
    //! Default Constructor.
    MeshContainer()
    { /* ... */ }

    //! Constructor.
    MeshContainer( const Teuchos::ArrayRCP<handle_type>& nodes,
		   const Teuchos::ArrayRCP<const double>& coords,
		   const int element_type,
		   const int element_topology,
		   const int nodes_per_element,
		   const Teuchos::ArrayRCP<handle_type>& elements,
		   const Teuchos::ArrayRCP<const handle_type>& connectivity )
	: d_nodes( nodes )
	, d_coords( coords )
	, d_element_type( element_type )
	, d_element_topology( element_topology )
	, d_nodes_per_element( nodes_per_element )
	, d_elements( elements )
	, d_connectivity( connectivity )
    { /* ... */ }

    //! Destructor.
    ~MeshContainer()
    { /* ... */ }

    //! Get the beginning of the nodes set.
    typename Teuchos::ArrayRCP<handle_type>::const_iterator nodesBegin() const
    { return d_nodes.begin(); }

    //! Get the end of the nodes set.
    typename Teuchos::ArrayRCP<handle_type>::const_iterator nodesEnd() const
    { return d_nodes.end(); }

    //! Get the beginning of the coordinates vector.
    Teuchos::ArrayRCP<const double>::const_iterator coordsBegin() const
    { return d_coords.begin(); }

    //! Get the end of the coordinates vector.
    Teuchos::ArrayRCP<const double>::const_iterator coordsEnd() const
    { return d_coords.end(); }

    //! Get the element type.
    std::size_t getElementType() const
    { return d_element_type; }

    //! Get the element topology.
    std::size_t getElementTopology() const
    { return d_element_topology; }

    //! Get the number of nodes constructing a single element.
    std::size_t getNodesPerElement() const
    { return d_nodes_per_element; }

    //! Get the beginning of the elements set.
    typename Teuchos::ArrayRCP<handle_type>::const_iterator elementsBegin() const
    { return d_elements.begin(); }

    //! Get the end of the elements set.
    typename Teuchos::ArrayRCP<handle_type>::const_iterator elementsEnd() const
    { return d_elements.end(); }

    //! Get the beginning of the connectivity vector.
    typename Teuchos::ArrayRCP<const handle_type>::const_iterator 
    connectivityBegin() const
    { return d_connectivity.begin(); }

    //! Get the ending of the connectivity vector.
    typename Teuchos::ArrayRCP<const handle_type>::const_iterator 
    connectivityEnd() const
    { return d_connectivity.end(); }
    
  private:

    // Nodes.
    Teuchos::ArrayRCP<handle_type> d_nodes;

    // Coordinates.
    Teuchos::ArrayRCP<const double> d_coords;

    // Element type.
    std::size_t d_element_type;

    // Element topology.
    std::size_t d_element_topology;

    // Nodes per element.
    std::size_t d_nodes_per_element;

    // Elements.
    Teuchos::ArrayRCP<handle_type> d_elements;

    // Connectivity.
    Teuchos::ArrayRCP<const handle_type> d_connectivity;
};

//---------------------------------------------------------------------------//
// MeshTraits specialization for the mesh container.
//---------------------------------------------------------------------------//
template<>
template<typename Handle>
struct MeshTraits< MeshContainer<Handle> >
{
    typedef MeshContainer<Handle> Container;
    typedef typename Container::handle_type handle_type;

    typedef typename Teuchos::ArrayRCP<Handle>::const_iterator 
    const_node_iterator;
    typedef Teuchos::ArrayRCP<const double>::const_iterator 
    const_coordinate_iterator;
    typedef typename Teuchos::ArrayRCP<Handle>::const_iterator 
    const_element_iterator;
    typedef typename Teuchos::ArrayRCP<const Handle>::const_iterator 
    const_connectivity_iterator;

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


    static inline std::size_t elementType( const Container& container )
    { return container.getElementType(); }

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
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit


#endif // end DTK_MESHCONTAINER_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshContainer.hpp
//---------------------------------------------------------------------------//
