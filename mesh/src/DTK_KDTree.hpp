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

#include "Teuchos_RCP.hpp"

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

    // Build the KDTree.
    void build();

    // Search the KDTree.
    ElementHandle findPoint( &const double coords[3] );
};

} // end namespace DataTransferKit

#endif // DTK_KDTREE_HPP

//---------------------------------------------------------------------------//
// end KDTree.hpp
//---------------------------------------------------------------------------//

