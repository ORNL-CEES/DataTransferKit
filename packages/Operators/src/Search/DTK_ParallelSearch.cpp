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
 * \brief DTK_ParallelSearch.cpp
 * \author Stuart R. Slattery
 * \brief Parallel search.
 */
//---------------------------------------------------------------------------//

#include "DTK_ParallelSearch.hpp"
#include "DTK_DBC.hpp"

#include <Tpetra_Distributor.hpp>

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
    : d_comm( comm )
    , d_physical_dim( physical_dimension )
    , d_track_missed_range_entities( false )
    , d_missed_range_entity_ids( 0 )
{
    // Determine if we are tracking missed range entities.
    if ( parameters.isParameter("Track Missed Range Entities") )
    {
	d_track_missed_range_entities =
	    parameters.get<bool>("Track Missed Range Entities");
    }

    // Build a coarse global search as this object must be collective across
    // the communicator.
    d_coarse_global_search = Teuchos::rcp(
	new CoarseGlobalSearch(d_comm, physical_dimension, 
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
// Search the domain with the range entity centroids and construct the
// graph. This will update the state of the object.
void ParallelSearch::search( 
    const EntityIterator& range_iterator,
    const Teuchos::RCP<EntityLocalMap>& range_local_map,
    const Teuchos::ParameterList& parameters )
{
    // Empty range flag.
    d_empty_range = ( 0 == range_iterator.size() );

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

    // If needed, extract the range entities that were missed during the
    // coarse global search.
    Teuchos::Array<EntityId> found_range_entity_ids;
    Teuchos::Array<int> found_range_ranks;
    Teuchos::Array<EntityId> missed_range_entity_ids;
    Teuchos::Array<int> missed_range_ranks;
    if ( d_track_missed_range_entities )
    {
	missed_range_entity_ids = Teuchos::Array<EntityId>( 
	    d_coarse_global_search->getMissedRangeEntityIds() );
	missed_range_ranks.assign( missed_range_entity_ids.size(),
				   d_comm->getRank() );
    }

    // Only do the local search if there are local domain entities.
    Teuchos::Array<int> export_range_ranks;
    Teuchos::Array<EntityId> export_data;
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
		// Store the range data in the domain parallel decomposition.
		std::pair<EntityId,int> range_rank(
		    range_entity_ids[n], range_owner_ranks[n] );

		std::pair<EntityId,EntityId> domain_range(
		    domain_parents[p].id(), range_entity_ids[n] );

		local_coords().assign( 
		    reference_coordinates(d_physical_dim*p,d_physical_dim) );
		std::pair<EntityId,Teuchos::Array<double> >
		    domain_ref_pair( domain_parents[p].id(), local_coords );

		d_range_owner_ranks.insert( range_rank );
		d_domain_to_range_map.insert( domain_range );
		ref_map.insert( domain_ref_pair );

		// Extract the data to communicate back to the range parallel
		// decomposition. 
		export_range_ranks.push_back( range_owner_ranks[n] );
		export_data.push_back( range_entity_ids[n] );
		export_data.push_back( domain_parents[p].id() );
		export_data.push_back( 
		    Teuchos::as<EntityId>(d_comm->getRank()) );
	    }

	    // If we found parents for the point, store them.
	    if ( num_parents > 0 )
	    {
		std::pair<EntityId,
			  std::unordered_map<EntityId,Teuchos::Array<double> > 
			  > range_ref_pair( range_entity_ids[n], ref_map );
		d_parametric_coords.insert( range_ref_pair );

		// If we are tracking missed entities, also track those that
		// we found so we can determine if an entity was found after
		// being sent to multiple destinations.
		if ( d_track_missed_range_entities )
		{
		    found_range_entity_ids.push_back( range_entity_ids[n] );
		    found_range_ranks.push_back( range_owner_ranks[n] );
		}
	    }
	    
	    // Otherwise, if we are tracking missed entities report this.
	    else if ( d_track_missed_range_entities )
	    {
		missed_range_entity_ids.push_back( range_entity_ids[n] );
		missed_range_ranks.push_back( range_owner_ranks[n] );
	    }
	}
    }

    // Back-communicate the domain entities in which we found each range
    // entity to complete the mapping.
    Tpetra::Distributor domain_to_range_dist( d_comm );
    int num_import = 
	domain_to_range_dist.createFromSends( export_range_ranks() );
    Teuchos::Array<EntityId> domain_data( 3.0*num_import );
    Teuchos::ArrayView<const EntityId> export_data_view = export_data();
    domain_to_range_dist.doPostsAndWaits( export_data_view, 3, domain_data() );

    // Store the domain data in the range parallel decomposition.
    for ( int i = 0; i < num_import; ++i )
    {
	std::pair<EntityId,int> domain_rank(
	    domain_data[3*i+1], domain_data[3*i+2] );

	std::pair<EntityId,EntityId> range_domain(
	    domain_data[3*i], domain_data[3*i+1] );

	d_domain_owner_ranks.insert( domain_rank );
	d_range_to_domain_map.insert( range_domain );
    }

    // If we are tracking missed entities, back-communicate the missing entities
    // and found entities to determine which entities are actually missing.
    if ( d_track_missed_range_entities )
    {
	// Back-communicate the missing entities.
	Tpetra::Distributor missed_range_dist( d_comm );
	int num_import_missed = 
	    missed_range_dist.createFromSends( missed_range_ranks() );
	Teuchos::Array<EntityId> import_missed( num_import_missed );
	Teuchos::ArrayView<const EntityId> missed_view = 
	    missed_range_entity_ids();
	missed_range_dist.doPostsAndWaits( missed_view, 1, import_missed() );

	// Back-communicate the found entities.
	Tpetra::Distributor found_range_dist( d_comm );
	int num_import_found = 
	    found_range_dist.createFromSends( found_range_ranks() );
	Teuchos::Array<EntityId> import_found( num_import_found );
	Teuchos::ArrayView<const EntityId> found_view = 
	    found_range_entity_ids();
	found_range_dist.doPostsAndWaits( found_view, 1, import_found() );

	// Intersect the found and missed entities to determine if there are any
	// that were found on one process but missed on another.
	std::sort( import_missed.begin(), import_missed.end() );
	std::sort( import_found.begin(), import_found.end() );
	Teuchos::Array<EntityId> false_positive_missed(
	    import_missed.size() + import_found.size() );
	auto false_positive_end = 
	    std::set_intersection( import_missed.begin(), import_missed.end(),
				   import_found.begin(), import_found.end(),
				   false_positive_missed.begin() );

	// Create a list of missed entities without the false positives.
	d_missed_range_entity_ids.resize( num_import_missed );
	auto missed_range_end = std::set_difference( 
	    import_missed.begin(), import_missed.end(),
	    false_positive_missed.begin(), false_positive_end,
	    d_missed_range_entity_ids.begin() );

	// Create a unique list of missed entities without the false positives.
	std::sort( d_missed_range_entity_ids.begin(), missed_range_end );
	auto missed_range_unique_end = std::unique(
	    d_missed_range_entity_ids.begin(), missed_range_end );
	d_missed_range_entity_ids.resize(
	    std::distance(d_missed_range_entity_ids.begin(),
			  missed_range_unique_end) );
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
    auto range_it = range_ids.begin();
    for ( auto domain_it = domain_pair.first;
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
    DTK_REQUIRE( !d_empty_range );
    auto range_pair = d_range_to_domain_map.equal_range( range_id );
    domain_ids.resize( std::distance(range_pair.first,range_pair.second) );
    auto domain_it = domain_ids.begin();
    for ( auto range_it = range_pair.first;
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
// Get the owner rank of a given domain entity.
int ParallelSearch::domainEntityOwnerRank( const EntityId domain_id ) const
{
    DTK_REQUIRE( !d_empty_range );
    DTK_REQUIRE( d_domain_owner_ranks.count(domain_id) );
    return d_domain_owner_ranks.find(domain_id)->second;
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
// Return the ids of the range entities that were not mapped during the last
// setup phase (i.e. those that are guaranteed to not receive data from the
// transfer). 
Teuchos::ArrayView<const EntityId> 
ParallelSearch::getMissedRangeEntityIds() const
{
    return d_missed_range_entity_ids();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_ParallelSearch.cpp
//---------------------------------------------------------------------------//
