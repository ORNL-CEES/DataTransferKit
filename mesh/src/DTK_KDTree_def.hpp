//---------------------------------------------------------------------------//
/*!
 * \file DTK_KDTree_def.hpp
 * \author Stuart R. Slattery
 * \brief KDTree definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_KDTREE_DEF_HPP
#define DTK_KDTREE_DEF_HPP

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

} // end namespace DataTransferKit

#endif // end DTK_KDTREE_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_KDTree_def.hpp
//---------------------------------------------------------------------------//

