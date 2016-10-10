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
 * \brief DTK_POD_PointCloudEntitySet.cpp
 * \author Stuart R. Slattery
 * \brief POD point cloud entity set.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "DTK_POD_PointCloudEntitySet.hpp"
#include "DTK_POD_PointCloudEntity.hpp"
#include "DTK_POD_PointCloudEntityIterator.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_DefaultMpiComm.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
POD_PointCloudEntitySet::POD_PointCloudEntitySet(
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const double* cloud_coords,
    const EntityId* global_ids,
    const unsigned num_points,
    const int space_dim,
    const DataLayout layout )
    : d_comm( comm )
    , d_cloud_coords( cloud_coords )
    , d_global_ids( global_ids )
    , d_num_points( num_points )
    , d_space_dim( space_dim )
    , d_layout( layout )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the parallel communicator for the entity set.
Teuchos::RCP<const Teuchos::Comm<int> >
POD_PointCloudEntitySet::communicator() const
{
    return d_comm;
}

//---------------------------------------------------------------------------//
// Return the largest physical dimension of the entities in the set. 
int POD_PointCloudEntitySet::physicalDimension() const
{
    return d_space_dim;
}

//---------------------------------------------------------------------------//
// Given an EntityId, get the entity.
void POD_PointCloudEntitySet::getEntity( const EntityId entity_id,
                                         const int topological_dimension,
                                         Entity& entity ) const
{
    DTK_REQUIRE( 0 == topological_dimension );
    DTK_REQUIRE( 1 == std::count(&d_global_ids[0],
                                 &d_global_ids[0] + d_num_points,
                                 entity_id) );
    
    int local_id = std::distance( &d_global_ids[0],
                                  std::find(&d_global_ids[0],
                                            &d_global_ids[0] + d_num_points,
                                            entity_id) );
    
    entity = POD_PointCloudEntity( d_cloud_coords,
                                   d_num_points,
                                   d_space_dim,
                                   d_layout,
                                   entity_id,
                                   local_id,
                                   d_comm->getRank() );
}

//---------------------------------------------------------------------------//
// Get an iterator over a subset of the entity set that satisfies the given
// predicate. 
EntityIterator POD_PointCloudEntitySet::entityIterator(
    const int topological_dimension,
    const PredicateFunction& predicate ) const
{
    EntityIterator iterator;
    
    if ( 0 == topological_dimension )
    {
        iterator = POD_PointCloudEntityIterator( d_cloud_coords,
                                                 d_global_ids,
                                                 d_num_points,
                                                 d_space_dim,
                                                 d_layout,
                                                 d_comm->getRank(),
                                                 predicate );
    }

    return iterator;
}

//---------------------------------------------------------------------------//
// Given an entity, get the entities of the given type that are adjacent to
// it. 
void POD_PointCloudEntitySet::getAdjacentEntities(
    const Entity& entity,
    const int adjacent_dimension,
    Teuchos::Array<Entity>& adjacent_entities ) const
{
    // There is no adjacency information in a point cloud.
    DTK_INSIST( false );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_POD_PointCloudEntitySet.cpp
//---------------------------------------------------------------------------//
