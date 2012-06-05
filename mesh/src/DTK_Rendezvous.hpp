//---------------------------------------------------------------------------//
/*!
 * \file DTK_Rendezvous.hpp
 * \author Stuart R. Slattery
 * \brief Rendezous declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_RENDEZVOUS_HPP
#define DTK_RENDEZVOUS_HPP

#include <vector>

#include "DTK_Mesh.hpp"
#include "DTK_KDTree.hpp"
#include "DTK_RCB.hpp"
#include "DTK_BoundingBox.hpp"
#include "DTK_Node.hpp"
#include "DTK_Element.hpp"
#include <DTK_NodeTraits.hpp>
#include <DTK_ElementTraits.hpp>
#include <DTK_FieldTraits.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

template<typename SourceMeshNodeField, typename SourceMeshElementField>
class Rendezvous
{
    
  public:

    //@{
    //! Typedefs.
    typedef SourceMeshNodeField                      source_node_field_type;
    typedef SourceMeshElementField                   source_element_field_type;
    typedef SourceMeshNodeField::value_type          source_node_type;
    typedef source_node_type::handle_type            source_node_handle_type;
    typedef source_node_type::coordinate_type        source_node_coordinate_type;
    typedef SourceMeshElementField::value_type       source_element_type;
    typedef source_element_type::handle_type         source_element_handle_type;
    typedef std::vector<source_element_handle_type>  SourceElementHandleVec;
    typedef Mesh<source_element_handle_type>         MeshType;
    typedef Teuchos::RCP<MeshType>                   RCP_Mesh;
    typedef KDTree<source_element_handle_type>       KDTreeType;
    typedef Teuchos::RCP<KDTreeType>                 RCP_KDTree;
    typedef RCB<SourceMeshNodeField>                 RCBType;
    typedef Teuchos::RCP<RCBType>                    RCP_RCB;
    typedef Teuchos::Comm<int>                       CommType;
    typedef Teuchos::RCP<const CommType>             RCP_Comm;
    //@}

    //@{
    //! Concrete node and element types.
    typedef Node<source_node_handle_type, source_node_coordinate_type> NodeType;
    typedef Element<
	source_element_handle_type,
	ElementTraits<source_element_type>::type(),
	ElementTraits<source_element_type>::topology(),
	ElementTraits<source_element_type>::numNodes() > ElementType;
    //@}

    // Constructor.
    Rendezvous( const RCP_Comm& global_comm, const BoundingBox& global_box );

    // Destructor.
    ~Rendezvous();

    // Build the rendezvous decomposition.
    void build( const SourceMeshNodeField& source_nodes,
		const SourceMeshElementField& source_elements );

    // Get the rendezvous process for a list of node coordinates.
    std::vector<int> 
    getRendezvousProc( const std::vector<double> &coords ) const;

    // Get the source mesh elements containing a list of coordinates.
    SourceElementHandleVec 
    getElementHandles( const std::vector<double>& coords ) const;

  private:

    // Extract the mesh nodes and elements that are in the bounding box.
    void getMeshInBox( const SourceMeshNodeField& source_nodes,
		       const SourceMeshElementField& source_elements,
		       std::vector<NodeType>& local_nodes,
		       std::vector<ElementType>& local_elements );

  private:

    // Global communicator over which to perform the rendezvous.
    RCP_Comm d_global_comm;

    // Bounding box over which to perform rendezvous.
    BoundingBox d_global_box;

    // RCB Partitioning.
    RCP_RCB d_rcb;

    // Local kD-tree.
    RCP_KDTree d_kdtree;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_Rendezvous_def.hpp"

#endif // end DTK_RENDEZVOUS_HPP

//---------------------------------------------------------------------------//
// end DTK_Rendezvous.hpp
//---------------------------------------------------------------------------//

