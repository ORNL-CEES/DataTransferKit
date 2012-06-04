//---------------------------------------------------------------------------//
/*!
 * \file DTK_Rendezvous.hpp
 * \author Stuart R. Slattery
 * \brief Rendezous declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_RENDEZVOUS_HPP
#define DTK_RENDEZVOUS_HPP

#include "DTK_Mesh.hpp"
#include "DTK_KDTree.hpp"
#include "DTK_RCB.hpp"
#include "DTK_BoundingBox.hpp"
#include <DTK_NodeTraits.hpp>
#include <DTK_ElementTraits.hpp>
#include <DTK_FieldTraits.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace DataTransferKit
{

template<typename SourceNodeField, typename SourceElementField,
	 typename TargetNodeField>
class Rendezvous
{
    
  public:

    //@{
    //! Typedefs.
    typedef SourceNodeField                     source_node_field_type;
    typedef SourceElementField                  source_element_field_type;
    typedef TargetNodeField                     target_node_field_type;
    typedef SourceNodeField::value_type         source_node_type;
    typedef SourceElementField::value_type      source_element_type;
    typedef source_element_type::handle_type    source_element_handle_type;
    typedef TargetNodeField::value_type         target_node_type;
    typedef Mesh<source_element_handle_type>    MeshType;
    typedef Teuchos::RCP<MeshType>              RCP_Mesh;
    typedef KDTree<source_element_handle_type>  KDTreeType;
    typedef Teuchos::RCP<KDTreeType>            RCP_KDTree;
    typedef RCB<SourceNodeField>                RCBType;
    typedef Teuchos::RCP<RCBType>               RCP_RCB;
    typedef Teuchos::Comm<int>                  CommType;
    typedef Teuchos::RCP<const CommType>        RCP_Comm;
    //@}

    // Constructor.
    Rendezvous( const RCP_Comm& global_comm );

    // Destructor.
    ~Rendezvous();

    // Build the rendezvous decomposition.
    void build( const SourceNodeField& source_nodes,
		const SourceElementField& source_elements,
		const TargetNodeField& target_nodes );

  private:

    // Global communicator over which to perform the rendezvous.
    RCP_Comm d_global_comm;

    // Local rendezvous mesh.
    RCP_Mesh d_mesh;

    // Local kD-tree.
    RCP_KDTree d_tree;
};

} // end namespace DataTransferKit

#endif // end DTK_RENDEZVOUS_HPP

//---------------------------------------------------------------------------//
// end DTK_Rendezvous.hpp
//---------------------------------------------------------------------------//

