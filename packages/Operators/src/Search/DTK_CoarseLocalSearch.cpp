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
 * \file DTK_CoarseLocalSearch.cpp
 * \author Stuart R. Slattery
 * \brief CoarseLocalSearch definition.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <algorithm>

#include "DTK_CoarseLocalSearch.hpp"
#include "DTK_DBC.hpp"
#include "DTK_SearchTreeFactory.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
CoarseLocalSearch::CoarseLocalSearch( 
    const EntityIterator& entity_iterator,
    const Teuchos::RCP<EntityLocalMap>& local_map,
    const Teuchos::ParameterList& parameters )
{
    // Setup the centroid array. These will be interleaved.
    int space_dim = 0;
    int num_entity = entity_iterator.size();
    if ( num_entity > 0 )
    {
	space_dim = entity_iterator.begin()->physicalDimension();
    }
    d_entity_centroids.resize( space_dim * num_entity );

    // Add the centroids.
    EntityIterator entity_it;
    EntityIterator begin_it = entity_iterator.begin();
    EntityIterator end_it = entity_iterator.end();
    int entity_local_id = 0;
    for ( entity_it = begin_it;
	  entity_it != end_it;
	  ++entity_it )
    {
	local_map->centroid( 
	    *entity_it, 
	    d_entity_centroids(space_dim*entity_local_id,space_dim) );
	d_entity_map.insert( 
	    std::pair<int,Entity>(entity_local_id,*entity_it) );
	++entity_local_id;
    }

    // Build a static search tree.
    int leaf_size = 20;
    if ( parameters.isParameter("Coarse Local Search Leaf Size") )
    {
	leaf_size = parameters.get<int>("Coarse Local Search Leaf Size");
    }
    leaf_size = std::min( leaf_size, num_entity );
    d_tree = SearchTreeFactory::createStaticTree(
	space_dim, d_entity_centroids(), leaf_size );
    DTK_ENSURE( Teuchos::nonnull(d_tree) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Find the set of entities a point neighbors.
 */
void CoarseLocalSearch::search( const Teuchos::ArrayView<const double>& point,
				const Teuchos::ParameterList& parameters,
				Teuchos::Array<Entity>& neighbors ) const
{
    // Find the leaf of nearest neighbors.
    int num_neighbors = 100;
    if ( parameters.isParameter("Coarse Local Search kNN") )
    {
	num_neighbors = parameters.get<int>("Coarse Local Search kNN");
    }
    num_neighbors = 
	std::min( num_neighbors, Teuchos::as<int>(d_entity_map.size()) );
    Teuchos::Array<unsigned> local_neighbors = 
	d_tree->nnSearch( point, num_neighbors );

    // Extract the neighbors.
    neighbors.resize( local_neighbors.size() );
    Teuchos::Array<unsigned>::const_iterator local_it;
    Teuchos::Array<Entity>::iterator entity_it;
    for ( local_it = local_neighbors.begin(),
	 entity_it = neighbors.begin();
	  local_it != local_neighbors.end();
	  ++local_it, ++entity_it )
    {
	DTK_CHECK( d_entity_map.count(*local_it) );
	*entity_it = d_entity_map.find(*local_it)->second;
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_CoarseLocalSearch.cpp
//---------------------------------------------------------------------------//

