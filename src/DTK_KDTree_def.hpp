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
#include "DTK_DBC.hpp"
#include "DataTransferKit_config.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_as.hpp>

#include <MBRange.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param mesh The RendezvousMesh over which to build the kD-tree. 
 *
 * \param dim The dimension of the kD-tree. This should be the same as the
 * dimension mesh.
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
    DTK_REMEMBER( moab::ErrorCode error );

    moab::Range elements;
#if HAVE_DTK_DBC
    error = d_mesh->getMoab()->get_entities_by_dimension( 0, d_dim, elements );
#else
    d_mesh->getMoab()->get_entities_by_dimension( 0, d_dim, elements );
#endif
    DTK_CHECK( moab::MB_SUCCESS == error );

#if HAVE_DTK_DBC
    error = d_tree.build_tree( elements , d_root );
#else
    d_tree.build_tree( elements , d_root );
#endif
    DTK_CHECK( moab::MB_SUCCESS == error );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find a point in the tree. Return false if we didn't find it in the
 * tree.
 *
 * \param coords Point coordinates to locate in the tree. Point dimensions
 * less than or equal to 3 are valid but the point most be the same dimension
 * as the tree.
 *
 * \param element The global ordinal of the client element the point was found
 * in. This global ordinal is not valid if this function returns false.
 *
 * \param tolerance Absolute tolerance for point searching. Will be used when
 * checking the reference cell ( and is therefore absolute ).
 *
 * \return Return true if the point was found in the kD-tree, false if not.
 */
template<typename GlobalOrdinal>
bool KDTree<GlobalOrdinal>::findPoint( const Teuchos::Array<double>& coords,
				       GlobalOrdinal& element,
				       double tolerance )
{
    DTK_REQUIRE( 0 <= coords.size() && coords.size() <= 3 );
    DTK_REQUIRE( (int) coords.size() == d_dim );

    double point[3];
    for ( int d = 0; d < d_dim; ++d )
    {
	point[d] = coords[d];
    }
    for ( int d = d_dim; d < 3; ++d )
    {
	point[d] = 0.0;
    }

    DTK_REMEMBER( moab::ErrorCode error );
    moab::EntityHandle leaf;
#if HAVE_DTK_DBC
    error = d_tree.leaf_containing_point( d_root, point, leaf );
#else
    d_tree.leaf_containing_point( d_root, point, leaf );
#endif
    DTK_CHECK( moab::MB_SUCCESS == error );

    moab::EntityHandle mb_element = 0;
    bool point_in_leaf = findPointInLeaf( coords, leaf, mb_element, tolerance );
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
 * \brief Get all of the elements in the leaf containing a point. 
 *
 * \param coords Point coordinates to locate in the tree. Point dimensions
 * less than or equal to 3 are valid but the point most be the same dimension
 * as the tree.
 *
 * \param elements The global ordinals of the local client elements in the
 * leaf that the point was found in. 
 */
template<typename GlobalOrdinal>
void KDTree<GlobalOrdinal>::findLeaf( const Teuchos::Array<double>& coords,
				      Teuchos::Array<GlobalOrdinal>& elements )
{
    DTK_REQUIRE( 0 <= coords.size() && coords.size() <= 3 );
    DTK_REQUIRE( Teuchos::as<int>(coords.size()) == d_dim );

    double point[3];
    for ( int d = 0; d < d_dim; ++d )
    {
	point[d] = coords[d];
    }
    for ( int d = d_dim; d < 3; ++d )
    {
	point[d] = 0.0;
    }

    DTK_REMEMBER( moab::ErrorCode error );
    moab::EntityHandle leaf;
#if HAVE_DTK_DBC
    error = d_tree.leaf_containing_point( d_root, point, leaf );
#else
    d_tree.leaf_containing_point( d_root, point, leaf );
#endif
    DTK_CHECK( moab::MB_SUCCESS == error );

    // Get the elements in the leaf.
    std::vector<moab::EntityHandle> leaf_elements;
#if HAVE_DTK_DBC
    error = d_mesh->getMoab()->get_entities_by_dimension( 
	leaf, d_dim, leaf_elements );
#else
    d_mesh->getMoab()->get_entities_by_dimension( 
	leaf, d_dim, leaf_elements );
#endif
    DTK_CHECK( moab::MB_SUCCESS == error );

    // Extract the client ordinals.
    elements.resize( leaf_elements.size() );
    typename Teuchos::Array<GlobalOrdinal>::iterator element_iterator;
    std::vector<moab::EntityHandle>::const_iterator leaf_iterator;
    for ( leaf_iterator = leaf_elements.begin(), 
       element_iterator = elements.begin();
	  leaf_iterator != leaf_elements.end();
	  ++leaf_iterator, ++element_iterator )
    {
	*element_iterator = d_mesh->getNativeOrdinal( *leaf_iterator );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find a point in a leaf by doing point-in-volume on all elements in
 * the leaf. Return if the point was not found.
 *
 * \param coords Point coordinates to locate in the tree. Point dimensions
 * less than or equal to 3 are valid but the point most be the same dimension
 * as the tree.
 *
 * \param leaf The subset of mesh elements to search for the point in.
 *
 * \param element The leaf element the point was found in. This element is not
 * valid if this function returns false.
 *
 * \return Return true if the point was found in the leaf, false if not.
 */
template<typename GlobalOrdinal>
bool KDTree<GlobalOrdinal>::findPointInLeaf( 
    const Teuchos::Array<double>& coords, 
    const moab::EntityHandle leaf,
    moab::EntityHandle& element,
    double tolerance )
{
    DTK_REQUIRE( 0 <= coords.size() && coords.size() <= 3 );
    DTK_REQUIRE( Teuchos::as<int>(coords.size()) == d_dim );

    // Get the elements in the leaf.
    DTK_REMEMBER( moab::ErrorCode error );
    std::vector<moab::EntityHandle> leaf_elements;
#if HAVE_DTK_DBC
    error = d_mesh->getMoab()->get_entities_by_dimension( 
	leaf, d_dim, leaf_elements );
#else
    d_mesh->getMoab()->get_entities_by_dimension( 
	leaf, d_dim, leaf_elements );
#endif
    DTK_CHECK( moab::MB_SUCCESS == error );

    // Search the leaf elements with the point.
    Teuchos::Array<double> point( coords );
    std::vector<moab::EntityHandle>::const_iterator leaf_iterator;
    for ( leaf_iterator = leaf_elements.begin();
	  leaf_iterator != leaf_elements.end();
	  ++leaf_iterator )
    {
	if ( TopologyTools::pointInElement( 
		 point, *leaf_iterator, d_mesh->getMoab(), tolerance ) )
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

