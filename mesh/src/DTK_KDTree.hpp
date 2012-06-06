//---------------------------------------------------------------------------//
/*!
 * \file DTK_KDTree.hpp
 * \author Stuart R. Slattery
 * \brief KDTree declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_KDTREE_HPP
#define DTK_KDTREE_HPP

#include "DTK_RendezvousMesh.hpp"

#include <Teuchos_RCP.hpp>

#include <MBAdaptiveKDTree.hpp>

namespace DataTransferKit
{

template<typename Handle>
class KDTree
{
  public:

    //@{
    //! Typedefs.
    typedef Handle                               handle_type;
    typedef RendezvousMesh<Handle>               RendezvousMeshType;
    typedef Teuchos::RCP<RendezvousMeshType>     RCP_RendezvousMesh;
    //@}

    // Constructor.
    KDTree( const RCP_RendezvousMesh& mesh );

    // Destructor.
    ~KDTree();

    // Build the kD-tree.
    void build();

    // Find a point in the tree.
    handle_type findPoint( double coords[3] );

  private:

    // Find a point in a leaf.
    moab::EntityHandle findPointInLeaf( double coords[3],
					const moab::EntityHandle leaf );

  private:

    // Moab Mesh.
    RCP_RendezvousMesh d_mesh;

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

