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
 * \file DTK_CoarseGlobalSearch.cpp
 * \author Stuart R. Slattery
 * \brief CoarseGlobalSearch definition.
 */
//---------------------------------------------------------------------------//

#include <limits>

#include "DTK_CoarseGlobalSearch.hpp"

#include <Teuchos_CommHelpers.hpp>

#include <Tpetra_Distributor.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
CoarseGlobalSearch::CoarseGlobalSearch( 
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const int physical_dimension,
    const EntityIterator& domain_iterator,
    const Teuchos::ParameterList& parameters )
    : d_comm( comm )
    , d_space_dim( physical_dimension )
{
    // Assemble the local domain bounding box.
    Teuchos::Tuple<double,6> domain_box;
    assembleBoundingBox( domain_iterator, domain_box );

    // Gather the bounding boxes from all domains.
    int comm_size = d_comm->getSize();
    Teuchos::Array<double> all_bounds( 6*comm_size );
    Teuchos::gatherAll<int,double>(
    	*d_comm, 6, domain_box.getRawPtr(), 
    	all_bounds.size(), all_bounds.getRawPtr() );

    // Extract the bounding boxes.
    d_domain_boxes.resize( comm_size );
    int id = 0;
    for ( int n = 0; n < comm_size; ++n )
    {
    	id = 6*n;
    	d_domain_boxes[n] = Teuchos::tuple( 
    	    all_bounds[id], all_bounds[id+1], all_bounds[id+2],
    	    all_bounds[id+3], all_bounds[id+4], all_bounds[id+5] );
    }
}

//---------------------------------------------------------------------------//
// Destructor.
CoarseGlobalSearch::~CoarseGlobalSearch()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Redistribute a set of range entity centroid coordinates with their owner
// ranks to the owning domain process.
void CoarseGlobalSearch::search( const EntityIterator& range_iterator,
				 const Teuchos::RCP<EntityLocalMap>& range_local_map,
				 const Teuchos::ParameterList& parameters,
				 Teuchos::Array<EntityId>& range_entity_ids,
				 Teuchos::Array<int>& range_owner_ranks,
				 Teuchos::Array<double>& range_centroids ) const
{
    // Assemble the local range bounding box.
    Teuchos::Tuple<double,6> range_box;
    assembleBoundingBox( range_iterator, range_box );
	
    // Find the domain boxes it intersects with.
    Teuchos::Array<int> neighbor_ranks;
    Teuchos::Array<Teuchos::Tuple<double,6> > neighbor_boxes;
    int num_domains = d_domain_boxes.size();
    for ( int n = 0; n < num_domains; ++n )
    {
	if ( boxesIntersect(range_box,d_domain_boxes[n]) )
	{
	    neighbor_ranks.push_back(n);
	    neighbor_boxes.push_back( d_domain_boxes[n] );
	}
    }

    // For each local range entity, find the neighbors we should send it to.
    int num_neighbors = neighbor_boxes.size();
    EntityIterator range_begin = range_iterator.begin();
    EntityIterator range_end = range_iterator.end();
    EntityIterator range_it;
    Teuchos::Array<EntityId> send_ids;
    Teuchos::Array<int> send_ranks;
    Teuchos::Array<Teuchos::Array<double> > send_centroids(d_space_dim);
    Teuchos::Array<double> centroid(d_space_dim);
    for ( range_it = range_begin; range_it != range_end; ++range_it )
    {
	// Get the centroid.
	range_local_map->centroid( *range_it, centroid() );

	// Check the neighbors.
	for ( int n = 0; n < num_neighbors; ++n )
	{
	    // If the centroid is in the box, add it to the send list.
	    if ( pointInBox(centroid(),neighbor_boxes[n]) )
	    {
		send_ids.push_back( range_it->id() );
		send_ranks.push_back( neighbor_ranks[n] );
		for ( int d = 0; d < d_space_dim; ++d )
		{
		    send_centroids[d].push_back( centroid[d] );
		}
	    }
	}
    }
    int num_send = send_ranks.size();
    Teuchos::Array<int> range_ranks( num_send, d_comm->getRank() );

    // Create a distributor.
    Tpetra::Distributor distributor(d_comm);
    int num_range_import = distributor.createFromSends( send_ranks() );

    // Redistribute the range entity ids.
    Teuchos::ArrayView<const EntityId> send_ids_view = send_ids();
    range_entity_ids.resize( num_range_import );
    distributor.doPostsAndWaits( send_ids_view, 1, range_entity_ids() );

    // Redistribute the range entity owner ranks.
    Teuchos::ArrayView<const int> range_ranks_view = range_ranks();
    range_owner_ranks.resize( num_range_import );
    distributor.doPostsAndWaits( range_ranks_view, 1, range_owner_ranks() );

    // Redistribute the range entity centroids.
    range_centroids.resize( d_space_dim*num_range_import );
    Teuchos::Array<double> temp_receives( num_range_import );
    for ( int d = 0; d < d_space_dim; ++d )
    {
	Teuchos::ArrayView<const double> send_centroids_view = 
	    send_centroids[d]();
	distributor.doPostsAndWaits( send_centroids_view, 1, temp_receives() );
	for ( int n = 0; n < num_range_import; ++n )
	{
	    range_centroids[ d_space_dim*n + d ] = temp_receives[n];
	}
    }
}

//---------------------------------------------------------------------------//
// Assemble the local bounding box around an iterator.
void CoarseGlobalSearch::assembleBoundingBox( 
    const EntityIterator& entity_iterator,
    Teuchos::Tuple<double,6>& bounding_box ) const
{
    double max = std::numeric_limits<double>::max();
    bounding_box = Teuchos::tuple( max, max, max, -max, -max, -max );
    Teuchos::Tuple<double,6> entity_bounds;
    EntityIterator entity_begin = entity_iterator.begin();
    EntityIterator entity_end = entity_iterator.end();
    EntityIterator entity_it;
    for ( entity_it = entity_begin; entity_it != entity_end; ++entity_it )
    {
	entity_it->boundingBox( entity_bounds );
	for ( int n = 0; n < 3; ++n )
	{
	    bounding_box[n] = 
		std::min( bounding_box[n], entity_bounds[n] );
	    bounding_box[n+3] = 
		std::max( bounding_box[n+3], entity_bounds[n+3] );
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end CoarseGlobalSearch.cpp
//---------------------------------------------------------------------------//

