//---------------------------------------------------------------------------//
/*!
 * \file DTK_KDTree_def.hpp
 * \author Stuart R. Slattery
 * \brief KDTree definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_KDTREE_DEF_HPP
#define DTK_KDTREE_DEF_HPP

#include <vector>

#include "DTK_TopologyTools.hpp"
#include <DTK_Exception.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename Handle>
KDTree<Handle>::KDTree( const RCP_RendezvousMesh& mesh )
: d_mesh( mesh )
, d_tree( d_mesh->getMoab().get() )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<typename Handle>
KDTree<Handle>::~KDTree()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Build the kD-tree.
 */
template<typename Handle>
void KDTree<Handle>::build()
{ 
    moab::ErrorCode error;
    error = d_tree.build_tree( d_mesh->getElements(), d_root );
    testInvariant( moab::MB_SUCCESS == error,
		   "Failed to construct kD-tree." );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find a point in the tree. Return the native handle of the element it
 * was found in. 
 */
template<typename Handle>
Handle KDTree<Handle>::findPoint( double coords[3] )
{
    moab::ErrorCode error;
    moab::EntityHandle leaf;
    error = d_tree.leaf_containing_point( d_root, coords, leaf );
    testInvariant( moab::MB_SUCCESS == error,
		   "Failed to search kD-tree." );

    moab::EntityHandle element = findPointInLeaf( coords, leaf );
    return d_mesh->getNativeHandle( element );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find a point in a leaf. Throw a MeshException if the point was not
 * found.
 */
template<typename Handle>
moab::EntityHandle KDTree<Handle>::findPointInLeaf( 
    double coords[3], const moab::EntityHandle leaf )
{
    moab::ErrorCode error;

    // Get the largest dimension of the elements in the leaf.
    int dim;
    int num_elements_by_dim = 0;
    for ( int d = 0; d < 4; ++d )
    {
	error = d_mesh->getMoab()->get_number_entities_by_dimension(
	    leaf, d, num_elements_by_dim );
	testInvariant( moab::MB_SUCCESS == error,
		       "Failed to get number of elements by dimension." );

	if ( num_elements_by_dim > 0 )
	{
	    dim = d;
	}
    }

    // Get the elements in the leaf.
    std::vector<moab::EntityHandle> leaf_elements;
    error = d_mesh->getMoab()->get_entities_by_dimension( 
	leaf, dim, leaf_elements );
    testInvariant( moab::MB_SUCCESS == error,
		   "Failed to get leaf elements" );

    // Search the leaf elements with the point.
    std::vector<moab::EntityHandle>::const_iterator leaf_iterator;
    for ( leaf_iterator = leaf_elements.begin();
	  leaf_iterator != leaf_elements.end();
	  ++leaf_iterator )
    {
	if ( TopologyTools::pointInElement( 
		 coords, *leaf_iterator, d_mesh->getMoab() ) )
	{
	    return *leaf_iterator;
	}
    }

    // We didn't find the point so we'll throw an exception.
    throw PointNotFound( "Point not found in kD-tree." );

    // Return a 0 handle if no element was found.
    return 0;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_KDTREE_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_KDTree_def.hpp
//---------------------------------------------------------------------------//

