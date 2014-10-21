//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \brief DTK_ParallelSearch.cpp
 * \author Stuart R. Slattery
 * \brief Parallel search.
 */
//---------------------------------------------------------------------------//

#include "DTK_ParallelSearch.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
ParallelSearch::ParallelSearch( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const int physical_dimension,
    const EntityIterator& domain_iterator,
    const Teuchos::RCP<EntityLocalMap>& domain_local_map,
    const Teuchos::ParameterList& parameters )
    : d_physical_dim( physical_dimension )
{
    d_coarse_global_search = Teuchos::rcp(
	new CoarseGlobalSearch(comm, physical_dimension, 
			       domain_iterator, parameters) );
    d_coarse_local_search = Teuchos::rcp(
	new CoarseLocalSearch(domain_iterator, domain_local_map, parameters) );
    d_fine_local_search = Teuchos::rcp(
	new FineLocalSearch(domain_local_map) );
}

//---------------------------------------------------------------------------//
// Destructor.
ParallelSearch::~ParallelSearch()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Search the domain with the range entity centroids and construct the
// graph. This will update the state of the object.
void ParallelSearch::search( 
    const EntityIterator& range_iterator,
    const Teuchos::RCP<EntityLocalMap>& range_local_map,
    const Teuchos::ParameterList& parameters )
{
    // Perform a coarse global search to redistribute the range entities.
    Teuchos::Array<EntityId> range_entity_ids;
    Teuchos::Array<int> range_owner_ranks;
    Teuchos::Array<double> range_centroids;
    d_coarse_global_search( 
	range_iterator, range_local_map, parameters,
	range_entity_ids, range_owner_ranks, range_centroids );

    // For each range centroid, perform a local search.
    int num_range = range_entity_ids.size();
    Teuchos::Array<EntityId> domain_neighbors;
    Teuchos::Array<EntityId> domain_parents;
    Teuchos::Array<double> reference_coordinates;
    Teuchos::Array<double> local_coords( d_physical_dim );
    int num_neighbors = 0;
    int num_parents = 0;
    for ( int n = 0; n < num_range; ++n )
    {
	// Perform a coarse local search to get the nearest domain entities to
	// the point.
	d_coarse_local_search( range_centroids(d_physical_dim*n,d_physical_dim),
			       parameters,
			       domain_neighbors );
	num_neighbors = domain_neighbors.size();
	
	// Perform a fine local search to get the entities the point maps to.
	d_fine_local_search( domain_neighbors,
			     range_centroids(d_physical_dim*n,d_physical_dim),
			     parameters,
			     domain_parents,
			     reference_coordinates );

	// Store the mapped point.
	num_parents = domain_parents.size();
	for ( int p = 0; p < num_parents; ++p )
	{
	    std::pair<EntityId,int> range_rank(
		range_entity_ids[n], range_owner_ranks[n] );

	    std::pair<EntityId,EntityId> domain_range(
		domain_parents[p].id(), range_entity_ids[n] );

	    std::pairt<EntityId,EntityId> range_domain(
		range_entity_ids[n], domain_parents[p].id() );

	    local_coords.assign( 
		reference_coordinates(d_physical_dim*p,d_physical_dim) );
	    std::pair<std::pair<EntityId,EntityId>,Teuchos::Array<double> > 
		ref_pair( domain_range, local_coords );

	    d_range_owner_ranks.insert( range_rank );
	    d_domain_to_range_map.insert( domain_range );
	    d_range_to_domain_map.insert( range_domain);
	    d_parametric_coordinates.insert( ref_pair );
    }
}

//---------------------------------------------------------------------------//
// Given a domain entity id, get the ids of the range entities that mapped to it.
void ParallelSearch::getRangeIdsFromDomain( 
    const EntityId& domain_id,
    Teuchos::Array<EntityId>& range_ids ) const
{
    std::pair<std::unordered_multimap<EntityId,EntityId>::const_iterator,
	      std::unordered_multimap<EntityId,EntityId>::const_iterator >
	domain_pair = d_domain_to_domain_map.equal_range( domain_id );
    std::unordered_multimap<EntityId,EntityId>::const_iterator domain_it;
    range_ids.resize( std::distance(domain_pair.first,domain_pair.second) );
    Teuchos::Array<EntityId>::iterator range_it;
    for ( domain_it = domain_pair.first, range_it = range_ids.begin();
	  domain_it != domain_pair.second;
	  ++domain_it, ++range_it )
    {
	*range_it = domain_it->second;
    }
}

//---------------------------------------------------------------------------//
// Given a range entity id, get the ids of the domain entities that it mapped to.
void ParallelSearch::getDomainIdsFromRange( 
    const EntityId& range_id,
    Teuchos::Array<EntityId>& domain_ids ) const
{
    std::pair<std::unordered_multimap<EntityId,EntityId>::const_iterator,
	      std::unordered_multimap<EntityId,EntityId>::const_iterator >
	range_pair = d_range_to_domain_map.equal_range( range_id );
    std::unordered_multimap<EntityId,EntityId>::const_iterator range_it;
    domain_ids.resize( std::distance(range_pair.first,range_pair.second) );
    Teuchos::Array<EntityId>::iterator domain_it;
    for ( range_it = range_pair.first, domain_it = domain_ids.begin();
	  range_it != range_pair.second;
	  ++range_it, ++domain_it )
    {
	*domain_it = range_it->second;
    }
}

//---------------------------------------------------------------------------//
// Get the parametric coordinates of the range entities in the domain entities.
void ParallelSearch::rangeParametricCoordinates( 
    const std::pair<EntityId,EntityId>& domain_range_pair,
    Teuchos::ArrayView<const double>& parametric_coords ) const
{
    DTK_REQUIRE( d_parametric_coords.count(domain_range_pair) );
    parametric_coords =	d_parametric_coords.find(domain_range_pair)->second();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_PARALLELSEARCH_HPP

//---------------------------------------------------------------------------//
// end DTK_ParallelSearch.hpp
//---------------------------------------------------------------------------//
