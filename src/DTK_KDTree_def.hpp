//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
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
#include "DTK_Assertion.hpp"

#include <Teuchos_ArrayView.hpp>

#include <MBRange.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename GlobalOrdinal>
KDTree<GlobalOrdinal>::KDTree( const RCP_RendezvousMesh& mesh, const int dim )
  : d_mesh( mesh )
  , d_dim( dim )
  , d_tree( d_mesh->getMoab().get() )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<typename GlobalOrdinal>
KDTree<GlobalOrdinal>::~KDTree()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Build the kD-tree.
 */
template<typename GlobalOrdinal>
void KDTree<GlobalOrdinal>::build()
{ 
    moab::ErrorCode error;

    moab::Range elements;
    error = d_mesh->getMoab()->get_entities_by_dimension( 0, d_dim, elements );
    testInvariant( moab::MB_SUCCESS == error );

    error = d_tree.build_tree( elements , d_root );
    testInvariant( moab::MB_SUCCESS == error );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find a point in the tree. Return false if we didn't find it in the
 * tree.
 */
template<typename GlobalOrdinal>
bool KDTree<GlobalOrdinal>::findPoint( const Teuchos::Array<double>& coords,
				       GlobalOrdinal& element )
{
    testPrecondition( 0 <= coords.size() && coords.size() <= 3 );
    testPrecondition( (int) coords.size() == d_dim );

    double point[3];
    for ( int d = 0; d < d_dim; ++d )
    {
	point[d] = coords[d];
    }
    for ( int d = d_dim; d < 3; ++d )
    {
	point[d] = 0.0;
    }

    moab::ErrorCode error;
    moab::EntityHandle leaf;
    error = d_tree.leaf_containing_point( d_root, point, leaf );
    testInvariant( moab::MB_SUCCESS == error );

    moab::EntityHandle mb_element;
    bool point_in_leaf = findPointInLeaf( coords, leaf, mb_element );
    if ( point_in_leaf )
    {
	element = d_mesh->getNativeOrdinal( mb_element );
	return true;
    }
    else
    {
	return false;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find a point in a leaf. Return if the point was not found.
 */
template<typename GlobalOrdinal>
bool KDTree<GlobalOrdinal>::findPointInLeaf( 
    const Teuchos::Array<double>& coords, 
    const moab::EntityHandle leaf,
    moab::EntityHandle& element )
{
    testPrecondition( 0 <= coords.size() && coords.size() <= 3 );
    testPrecondition( (int) coords.size() == d_dim );

    moab::ErrorCode error;

    // Get the elements in the leaf.
    std::vector<moab::EntityHandle> leaf_elements;
    error = d_mesh->getMoab()->get_entities_by_dimension( 
	leaf, d_dim, leaf_elements );
    testInvariant( moab::MB_SUCCESS == error );

    // Search the leaf elements with the point.
    Teuchos::Array<double> point( coords );
    std::vector<moab::EntityHandle>::const_iterator leaf_iterator;
    for ( leaf_iterator = leaf_elements.begin();
	  leaf_iterator != leaf_elements.end();
	  ++leaf_iterator )
    {
	if ( TopologyTools::pointInElement( 
		 point, *leaf_iterator, d_mesh->getMoab() ) )
	{
	    element = *leaf_iterator;
	    return true;
	}
    }

    // Return false if no element was found.
    return false;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_KDTREE_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_KDTree_def.hpp
//---------------------------------------------------------------------------//

