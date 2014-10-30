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

    // Only do the local search if there are local domain entities.
    d_empty_domain = ( 0 == domain_iterator.size() );
    if ( !d_empty_domain )
    {
	d_coarse_local_search = Teuchos::rcp(
	    new CoarseLocalSearch(domain_iterator, domain_local_map, parameters) );
	d_fine_local_search = Teuchos::rcp(
	    new FineLocalSearch(domain_local_map) );
    }
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
    // Reset the state of the object.
    d_range_owner_ranks.clear();
    d_domain_to_range_map.clear();
    d_range_to_domain_map.clear();
    d_parametric_coords.clear();

    // Perform a coarse global search to redistribute the range entities.
    Teuchos::Array<EntityId> range_entity_ids;
    Teuchos::Array<int> range_owner_ranks;
    Teuchos::Array<double> range_centroids;
    d_coarse_global_search->search( 
	range_iterator, range_local_map, parameters,
	range_entity_ids, range_owner_ranks, range_centroids );

    // Only do the local search if there are local domain entities.
    if ( !d_empty_domain )
    {
	// For each range centroid, perform a local search.
	int num_range = range_entity_ids.size();
	Teuchos::Array<Entity> domain_neighbors;
	Teuchos::Array<Entity> domain_parents;
	Teuchos::Array<double> reference_coordinates;
	Teuchos::Array<double> local_coords( d_physical_dim );
	int num_parents = 0;
	for ( int n = 0; n < num_range; ++n )
	{
	    // Perform a coarse local search to get the nearest domain
	    // entities to the point.
	    d_coarse_local_search->search( 
		range_centroids(d_physical_dim*n,d_physical_dim),
		parameters,
		domain_neighbors );
	
	    // Perform a fine local search to get the entities the point maps
	    // to.
	    d_fine_local_search->search( 
		domain_neighbors,
		range_centroids(d_physical_dim*n,d_physical_dim),
		parameters,
		domain_parents,
		reference_coordinates );

	    // Store the potentially multiple parametric realizations of the
	    // point.
	    std::unordered_map<EntityId,Teuchos::Array<double> > ref_map;
	    num_parents = domain_parents.size();
	    for ( int p = 0; p < num_parents; ++p )
	    {
		std::pair<EntityId,int> range_rank(
		    range_entity_ids[n], range_owner_ranks[n] );

		std::pair<EntityId,EntityId> domain_range(
		    domain_parents[p].id(), range_entity_ids[n] );

		std::pair<EntityId,EntityId> range_domain(
		    range_entity_ids[n], domain_parents[p].id() );

		local_coords().assign( 
		    reference_coordinates(d_physical_dim*p,d_physical_dim) );
		std::pair<EntityId,Teuchos::Array<double> >
		    domain_ref_pair( domain_parents[p].id(), local_coords );

		d_range_owner_ranks.insert( range_rank );
		d_domain_to_range_map.insert( domain_range );
		d_range_to_domain_map.insert( range_domain );
		ref_map.insert( domain_ref_pair );
	    }
	    if ( num_parents > 0 )
	    {
		std::pair<EntityId,
			  std::unordered_map<EntityId,Teuchos::Array<double> > 
			  > range_ref_pair( range_entity_ids[n], ref_map );
		d_parametric_coords.insert( range_ref_pair );
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// Given a domain entity id, get the ids of the range entities that mapped to it.
void ParallelSearch::getRangeEntitiesFromDomain( 
    const EntityId domain_id,
    Teuchos::Array<EntityId>& range_ids ) const
{
    DTK_REQUIRE( !d_empty_domain );
    auto domain_pair = d_domain_to_range_map.equal_range( domain_id );
    range_ids.resize( std::distance(domain_pair.first,domain_pair.second) );
    for ( auto domain_it = domain_pair.first, auto range_it = range_ids.begin();
	  domain_it != domain_pair.second;
	  ++domain_it, ++range_it )
    {
	*range_it = domain_it->second;
    }
}

//---------------------------------------------------------------------------//
// Given a range entity id, get the ids of the domain entities that it mapped to.
void ParallelSearch::getDomainEntitiesFromRange( 
    const EntityId range_id,
    Teuchos::Array<EntityId>& domain_ids ) const
{
    DTK_REQUIRE( !d_empty_domain );
    auto range_pair = d_range_to_domain_map.equal_range( range_id );
    domain_ids.resize( std::distance(range_pair.first,range_pair.second) );
    for ( auto range_it = range_pair.first, auto domain_it = domain_ids.begin();
	  range_it != range_pair.second;
	  ++range_it, ++domain_it )
    {
	*domain_it = range_it->second;
    }
}

//---------------------------------------------------------------------------//
// Get the owner rank of a given range entity.
int ParallelSearch::rangeEntityOwnerRank( const EntityId range_id ) const
{
    DTK_REQUIRE( !d_empty_domain );
    DTK_REQUIRE( d_range_owner_ranks.count(range_id) );
    return d_range_owner_ranks.find(range_id)->second;
}

//---------------------------------------------------------------------------//
// Get the parametric coordinates of the range entities in the domain entities.
void ParallelSearch::rangeParametricCoordinatesInDomain( 
    const EntityId domain_id,
    const EntityId range_id,
    Teuchos::ArrayView<const double>& parametric_coords ) const
{
    DTK_REQUIRE( !d_empty_domain );
    DTK_REQUIRE( d_parametric_coords.count(range_id) );
    DTK_REQUIRE( d_parametric_coords.find(range_id)->second.count(domain_id) );
    parametric_coords =	
	d_parametric_coords.find(range_id)->second.find(domain_id)->second;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_ParallelSearch.cpp
//---------------------------------------------------------------------------//
