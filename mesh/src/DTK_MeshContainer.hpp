//---------------------------------------------------------------------------//
/*!
 * \file DTK_MeshContainer.hpp
 * \author Stuart R. Slattery
 * \brief A simple mesh container for rebuilding mesh data after serialization.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHCONTAINER_HPP
#define DTK_MESHCONTAINER_HPP

#include <iterator>
#include <vector>
#include <set>

#include <DTK_MeshTraits.hpp>

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
    MeshContainer( const std::set<handle_type>& nodes,
		   const std::vector<double>& coords,
		   const int element_type,
		   const int element_topology,
		   const int nodes_per_element,
		   const std::set<handle_type>& elements,
		   const std::vector<handle_type>& connectivity )
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
    typename std::set<handle_type>::const_iterator nodesBegin()
    { return d_nodes.begin(); }

    //! Get the end of the nodes set.
    typename std::set<handle_type>::const_iterator nodesEnd()
    { return d_nodes.end(); }

    //! Get the beginning of the coordinates vector.
    std::vector<double>::const_iterator coordsBegin()
    { return d_coords.begin(); }

    //! Get the end of the coordinates vector.
    std::vector<double>::const_iterator coordsEnd()
    { return d_coords.end(); }

    //! Get the element type.
    std::size_t getElementType()
    { return d_element_type; }

    //! Get the element topology.
    std::size_t getElementTopology()
    { return d_element_topology; }

    //! Get the number of nodes constructing a single element.
    std::size_t getNodesPerElement()
    { return d_nodes_per_element; }

    //! Get the beginning of the elements set.
    typename std::set<handle_type>::const_iterator elementsBegin()
    { return d_elements.begin(); }

    //! Get the end of the elements set.
    typename std::set<handle_type>::const_iterator elementsEnd()
    { return d_elements.end(); }

    //! Get the beginning of the connectivity vector.
    typename std::vector<handle_type>::const_iterator connectivityBegin()
    { return d_connectivity.begin(); }

    //! Get the endning of the connectivity vector.
    typename std::vector<handle_type>::const_iterator connectivityEnd()
    { return d_connectivity.end(); }
    
  private:

    // Nodes.
    std::set<handle_type> d_nodes;

    // Coordinates.
    std::vector<double> d_coords;

    // Element type.
    std::size_t d_element_type;

    // Element topology.
    std::size_t d_element_topology;

    // Nodes per element.
    std::size_t d_nodes_per_element;

    // Elements.
    std::set<handle_type> d_elements;

    // Connectivity.
    std::vector<handle_type> d_connectivity;
};

//---------------------------------------------------------------------------//
// MeshTraits specialization for the mesh container.
//---------------------------------------------------------------------------//
template<>
template<typename Handle>
struct MeshTraits< MeshContainer<Handle> >
{
    typedef MeshContainer<Handle> Container;

    static inline typename std::set<Handle>::const_iterator 
    nodesBegin( const Container& container )
    { return container.nodesBegin(); }

    static inline typename std::set<Handle>::const_iterator 
    nodesEnd( const Container& container )
    { return container.nodesEnd(); }

    static inline bool interleavedCoordinates( const Container& container )
    { return true; }

    static inline std::vector<double>::const_iterator 
    coordsBegin( const Container& container )
    { return container.coordsBegin(); }

    static inline std::vector<double>::const_iterator 
    coordsEnd( const Container& container )
    { return container.coordsEnd(); }


    static inline std::size_t elementType( const Container& container )
    { return container.getElementType(); }

    static inline std::size_t elementTopology( const Container& container )
    { return container.getElementTopology(); }

    static inline std::size_t nodesPerElement( const Container& container )
    { return container.getNodesPerElement(); }

    static inline typename std::set<Handle>::const_iterator 
    elementsBegin( const Container& container )
    { return container.elementsBegin(); }

    static inline typename std::set<Handle>::const_iterator 
    elementsEnd( const Container& container )
    { return container.elementsEnd(); }

    static inline typename std::vector<Handle>::const_iterator 
    connectivityBegin( const Container& container )
    { return container.connectivityBegin(); }

    static inline typename std::vector<Handle>::const_iterator 
    connectivityEnd( const Container& container )
    { return container.connectivityEnd(); }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit


#endif // end DTK_MESHCONTAINER_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshContainer.hpp
//---------------------------------------------------------------------------//
