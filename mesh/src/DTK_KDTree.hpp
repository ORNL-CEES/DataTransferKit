//---------------------------------------------------------------------------//
/*!
 * \file DTK_KDTree.hpp
 * \author Stuart R. Slattery
 * \brief KDTree declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_KDTREE_HPP
#define DTK_KDTREE_HPP

#include "DTK_Mesh.hpp"

#include <Teuchos_RCP.hpp>

#include <MBAdaptiveKDTree.hpp>

namespace DataTransferKit
{

template<typename ElementHandle>
class KDTree
{
  public:

    //@{
    //! Typedefs.
    typedef ElementHandle                      element_handle_type;
    typedef Mesh<ElementHandle>                MeshType;
    typedef Teuchos::RCP<MeshType>             RCP_Mesh;
    //@}

    // Constructor.
    KDTree( const RCP_Mesh& mesh );

    // Destructor.
    ~KDTree();

    // Build the kD-tree.
    void build();

    // Find a point in the tree.
    ElementHandle findPoint( double coords[3] );

  private:

    // Find a point in a leaf.
    moab::EntityHandle findPointInLeaf( double coords[3],
					const moab::EntityHandle leaf );

  private:

    // Mesh.
    RCP_Mesh d_mesh;

    // Adaptive kD-tree.
    moab::AdaptiveKDTree d_tree;

    // Tree root.
    moab::EntityHandle d_root;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_KDTree_def.hpp"

#endif // DTK_KDTREE_HPP

//---------------------------------------------------------------------------//
// end KDTree.hpp
//---------------------------------------------------------------------------//

