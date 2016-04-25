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
 * \file   DTK_SplineInterpolationPairing_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Local child/parent center pairings.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_SPLINEINTERPOLATIONPAIRING_IMPL_HPP
#define DTK_SPLINEINTERPOLATIONPAIRING_IMPL_HPP

#include "DTK_SplineInterpolationPairing.hpp"
#include "DTK_DBC.hpp"
#include "DTK_StaticSearchTree.hpp"
#include "DTK_EuclideanDistance.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<int DIM>
SplineInterpolationPairing<DIM>::SplineInterpolationPairing( 
    const Teuchos::ArrayView<const double>& child_centers,
    const Teuchos::ArrayView<const double>& parent_centers,
    const bool use_knn,
    const unsigned num_neighbors,
    const double radius )
{
    DTK_REQUIRE( 0 == child_centers.size() % DIM );
    DTK_REQUIRE( 0 == parent_centers.size() % DIM );

    // Setup the kD-tree
    unsigned leaf_size = 30;
    NanoflannTree<DIM> tree( child_centers, leaf_size );

    // Allocate arrays
    unsigned num_parents = parent_centers.size() / DIM;
    d_pairings.resize( num_parents );
    d_pair_sizes = Teuchos::ArrayRCP<EntityId>( num_parents );
    d_radii.resize( num_parents );

    // Search for pairs
    for ( unsigned i = 0; i < num_parents; ++i )
    {
	// if kNN do the nearest neighbor search for kNN+1 and calculate a
	// radius
	if ( use_knn )
	{
	    // Get the knn neighbors
	    d_pairings[i] =
		tree.nnSearch( parent_centers(DIM*i,DIM), num_neighbors + 1 );

	    // Get the radius from kNN + 1
	    d_radii[i] =
		EuclideanDistance<DIM>::distance(
		    parent_centers(DIM*i,DIM).getRawPtr(),
		    child_centers(DIM*d_pairings[i].back(),DIM).getRawPtr() );

	    // Pop the back of the neighbors to get kNN
	    if ( d_pairings[i].size() > 0 )
	    {
		d_pairings[i].pop_back();
	    }
	}

	// otherwise do the radius search
	else
	{
	    d_pairings[i] =
		tree.radiusSearch( parent_centers(DIM*i,DIM), radius );
	    d_radii[i] = radius;
	}

	// get the size of the support.
	d_pair_sizes[i] = d_pairings[i].size();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Given a parent center local id get the ids of the child centers
 * within the given radius.
 */
template<int DIM>
Teuchos::ArrayView<const unsigned> 
SplineInterpolationPairing<DIM>::childCenterIds(
    const unsigned parent_id ) const
{
    DTK_REQUIRE( parent_id < d_pairings.size() );
    Teuchos::ArrayView<const unsigned> id_view = d_pairings[parent_id]();
    return id_view;
}

//---------------------------------------------------------------------------//
/*
 * \brief Get the support radius of a given parent.
 */
template<int DIM>
double SplineInterpolationPairing<DIM>::parentSupportRadius(
    const unsigned parent_id ) const
{
    REQUIRE( parent_id < d_radii.size() );
    return d_radii[i];
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_SPLINEINTERPOLATIONPAIRING_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_SplineInterpolationPairing_impl.hpp
//---------------------------------------------------------------------------//

