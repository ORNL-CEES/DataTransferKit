//---------------------------------------------------------------------------//
/*!
 * \file DTK_KDTree_def.hpp
 * \author Stuart R. Slattery
 * \brief KDTree definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_KDTREE_DEF_HPP
#define DTK_KDTREE_DEF_HPP

#include <DTK_Exception.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename ElementHandle>
KDTree<ElementHandle>::KDTree( const RCP_Mesh& mesh )
: d_mesh( mesh )
, d_tree( mesh->getMoab()->get() )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<typename ElementHandle>
KDTree<ElementHandle>::~KDTree()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Build the tree.
 */
template<typename ElementHandle>
void KDTree<ElementHandle>::build()
{ 
    moab::ErrorCode error;
    error = d_tree.build_tree( mesh->getElements(), d_root );
    testInvariant( moab::MB_SUCCESS == error,
		   "Failed to construct kD-tree." );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find a point in the tree. Return the native handle of the element it
 * was found in. Throw a MeshException if the point was not found.
 */
template<typename ElementHandle>
ElementHandle KDTree<ElementHandle>::findPoint( const double coords[3] )
{
    moab::ErrorCode error;
    moab::EntityHandle leaf;
    error = d_tree.leaf_containing_point( d_root, coords, leaf );
    testInvariant( moab::MB_SUCCESS == error,
		   "Failed to search kD-tree." );

    moab::EntityHandle element;
    bool found_element = findPointInLeaf( coords, leaf, element );
    if ( !found element)
    {
	throw PointNotFound();
    }

    return d_mesh->getNativeHandle( element );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find a point in a leaf.
 */
bool KDTree<ElementHandle>::findPointInLeaf( const double coords[3],
					     const moab::EntityHandle leaf,
					     moab::EntityHandle &element )
{
    bool found_element = false;


}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_KDTREE_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_KDTree_def.hpp
//---------------------------------------------------------------------------//

