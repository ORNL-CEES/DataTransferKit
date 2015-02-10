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
 * \brief DTK_MoabEntityImpl.cpp
 * \author Stuart R. Slattery
 * \brief Moab entity implementation.
 */
//---------------------------------------------------------------------------//

#include <limits>

#include "DTK_MoabEntityImpl.hpp"
#include "DTK_MoabHelpers.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_Array.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
MoabEntityImpl::MoabEntityImpl(
    const moab::EntityHandle& moab_entity,
    const Teuchos::Ptr<moab::ParallelComm>& moab_mesh,
    const Teuchos::Ptr<MoabMeshSetIndexer>& set_indexer )
    : d_extra_data( new MoabEntityExtraData(moab_entity) )
    , d_moab_mesh( moab_mesh )
    , d_set_indexer( set_indexer )
{ /* ... */ }

//---------------------------------------------------------------------------//
//brief Destructor.
MoabEntityImpl::~MoabEntityImpl()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the entity type.
EntityType MoabEntityImpl::entityType() const
{
    return MoabHelpers::getEntityTypeFromMoabType(
	d_moab_mesh->get_moab()->type_from_handle(
	    d_extra_data->d_moab_entity) );
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
EntityId MoabEntityImpl::id() const
{ 
    return d_extra_data->d_moab_entity;
}
    
//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
int MoabEntityImpl::ownerRank() const
{ 
    int owner_rank = -1;
    DTK_CHECK_ERROR_CODE(
	d_moab_mesh->get_owner( d_extra_data->d_moab_entity, owner_rank )
	);
    return owner_rank;
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
int MoabEntityImpl::physicalDimension() const
{ 
    int dimension = 0;
    DTK_CHECK_ERROR_CODE(
	d_moab_mesh->get_moab()->get_dimension( dimension )
	);
    return dimension;
}

//---------------------------------------------------------------------------//
// Return the Cartesian bounding box around an entity.
void MoabEntityImpl::boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    Teuchos::Array<double> coordinates;

    // Node case.
    if ( ENTITY_TYPE_NODE == this->entityType() )
    {
	coordinates.resize( 3 );
	DTK_CHECK_ERROR_CODE(
	    d_moab_mesh->get_moab()->get_coords( 
		&(d_extra_data->d_moab_entity),
		1,
		coordinates.getRawPtr() )
	    );
    }

    // Element/face/edge case.
    else
    {
	MoabHelpers::getEntityNodeCoordinates( d_extra_data->d_moab_entity,
					       d_moab_mesh,
					       coordinates );
    }

    int num_nodes = coordinates.size() / 3;
    int space_dim = this->physicalDimension();
    double max = std::numeric_limits<double>::max();
    bounds = Teuchos::tuple( max, max, max, -max, -max, -max );
    for ( int d = 0; d < space_dim; ++d )
    {
	for ( int n = 0; n < num_nodes; ++n )
	{
	    bounds[d] = std::min( bounds[d], coordinates[n*space_dim + d] );
	    bounds[d+3] = std::max( bounds[d+3], coordinates[n*space_dim + d] );
	}
    }
}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
bool MoabEntityImpl::inBlock( const int block_id ) const
{
    moab::EntityHandle block_set = 
	d_set_indexer->getMeshSetFromIndex( block_id );
    return d_moab_mesh->get_moab()->contains_entities(
	block_set, &d_extra_data->d_moab_entity, 1 );
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
bool MoabEntityImpl::onBoundary( const int boundary_id ) const
{
    return this->inBlock( boundary_id );
}

//---------------------------------------------------------------------------//
// Get the extra data on the entity.
Teuchos::RCP<EntityExtraData> MoabEntityImpl::extraData() const
{
    return d_extra_data;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MoabEntityImpl.cpp
//---------------------------------------------------------------------------//
