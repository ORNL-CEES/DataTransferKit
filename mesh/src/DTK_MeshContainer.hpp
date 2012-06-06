//---------------------------------------------------------------------------//
/*!
 * \file DTK_MeshContainer.hpp
 * \author Stuart R. Slattery
 * \brief A simple mesh container for serializing and rebuilding mesh data.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MESHCONTAINER_HPP
#define DTK_MESHCONTAINER_HPP

#include <iterator>
#include <vector>
#include <set>

#include <DTK_MeshTraits.hpp>

#include <Teuchos_SerializationTraits.hpp>

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
    typedef MeshContainer<Handle> MeshType;

    static inline typename std::set<Handle>::const_iterator 
    nodesBegin( const MeshType& mesh )
    { return mesh.nodesBegin(); }

    static inline typename std::set<Handle>::const_iterator 
    nodesEnd( const MeshType& mesh )
    { return mesh.nodesEnd(); }

    static inline bool interleavedCoordinates( const MeshType& mesh )
    { return true; }

    static inline std::vector<double>::const_iterator 
    coordsBegin( const MeshType& mesh )
    { return mesh.coordsBegin(); }

    static inline std::vector<double>::const_iterator 
    coordsEnd( const MeshType& mesh )
    { return mesh.coordsEnd(); }


    static inline std::size_t elementType( const MeshType& mesh )
    { return mesh.getElementType(); }

    static inline std::size_t elementTopology( const MeshType& mesh )
    { return mesh.getElementTopology(); }

    static inline std::size_t nodesPerElement( const MeshType& mesh )
    { return mesh.getNodesPerElement(); }

    static inline typename std::set<Handle>::const_iterator 
    elementsBegin( const MeshType& mesh )
    { return mesh.elementsBegin(); }

    static inline typename std::set<Handle>::const_iterator 
    elementsEnd( const MeshType& mesh )
    { return mesh.elementsEnd(); }

    static inline typename std::vector<Handle>::const_iterator 
    connectivityBegin( const MeshType& mesh )
    { return mesh.connectivityBegin(); }

    static inline typename std::vector<Handle>::const_iterator 
    connectivityEnd( const MeshType& mesh )
    { return mesh.connectivityEnd(); }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Teuchos serialization traits specializations.
//---------------------------------------------------------------------------//
namespace Teuchos
{

template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::MeshContainer<int> >
    : public DirectSerializationTraits< Ordinal,
					DataTransferKit::MeshContainer<int> >
{ /* ... */ };

} // end namespace Teuchos

//---------------------------------------------------------------------------//

#endif // end DTK_MESHCONTAINER_HPP

//---------------------------------------------------------------------------//
// end DTK_MeshContainer.hpp
//---------------------------------------------------------------------------//
