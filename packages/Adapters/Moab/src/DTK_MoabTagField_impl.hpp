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
 * \brief DTK_MoabTagField_impl.hpp
 * \author Stuart R. Slattery
 * \brief Moab tag vector manager
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOABTAGFIELD_IMPL_HPP
#define DTK_MOABTAGFIELD_IMPL_HPP

#include <vector>

#include "DTK_DBC.hpp"
#include "DTK_MoabHelpers.hpp"

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template <class Scalar>
MoabTagField<Scalar>::MoabTagField(
    const Teuchos::RCP<moab::ParallelComm> &moab_mesh,
    const Teuchos::RCP<MoabMeshSetIndexer> set_indexer,
    const moab::EntityHandle &mesh_set, const moab::Tag &tag )
    : d_moab_mesh( moab_mesh )
    , d_set_indexer( set_indexer )
    , d_mesh_set( mesh_set )
    , d_tag( tag )
    , d_entity_dim( -1 )
{
    // Get the dimension of the tag.
    DTK_CHECK_ERROR_CODE(
        d_moab_mesh->get_moab()->tag_get_length( d_tag, d_tag_dim ) );

    // Get the entities in the set.
    std::vector<moab::EntityHandle> entities;
    DTK_CHECK_ERROR_CODE( d_moab_mesh->get_moab()->get_entities_by_handle(
        d_mesh_set, entities ); );

    // Get the global ids.
    int num_entities = entities.size();
    if ( 0 < num_entities )
    {
        // Get the topological dimension of the entities supporting the
        // field. We assume all entities supporting a field are of the same
        // type.
        d_entity_dim = MoabHelpers::getTopologicalDimensionFromMoabType(
            d_moab_mesh->get_moab()->type_from_handle( entities[0] ) );

        // Get the global ids of the entities supporting the field.
        Teuchos::Array<EntityId> global_ids( num_entities );
        MoabHelpers::getGlobalIds( *d_moab_mesh, entities.data(), num_entities,
                                   global_ids.getRawPtr() );

        // Create locally-owned support ids.
        int owner_rank = -1;
        int rank = d_moab_mesh->rank();
        for ( int n = 0; n < num_entities; ++n )
        {
            DTK_CHECK_ERROR_CODE(
                d_moab_mesh->get_owner( entities[n], owner_rank ) );
            if ( rank == owner_rank )
            {
                d_support_ids.push_back( global_ids[n] );
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Get the dimension of the field.
template <class Scalar>
int MoabTagField<Scalar>::dimension() const
{
    return d_tag_dim;
}

//---------------------------------------------------------------------------//
// Get the locally-owned support location ids of the field.
template <class Scalar>
Teuchos::ArrayView<const SupportId>
MoabTagField<Scalar>::getLocalSupportIds() const
{
    return d_support_ids();
}

//---------------------------------------------------------------------------//
// Given a local support id and a dimension, read data from the application
// field.
template <class Scalar>
double MoabTagField<Scalar>::readFieldData( const SupportId support_id,
                                            const int dimension ) const
{
    moab::EntityHandle entity =
        d_set_indexer->getEntityFromGlobalId( support_id, d_entity_dim );
    const void *tag_data = 0;
    DTK_CHECK_ERROR_CODE( d_moab_mesh->get_moab()->tag_get_by_ptr(
        d_tag, &entity, 1, &tag_data ) );
    return static_cast<const Scalar *>( tag_data )[dimension];
}

//---------------------------------------------------------------------------//
// Given a local support id, dimension, and field value, write data into the
// application field.
template <class Scalar>
void MoabTagField<Scalar>::writeFieldData( const SupportId support_id,
                                           const int dimension,
                                           const double data )
{
    moab::EntityHandle entity =
        d_set_indexer->getEntityFromGlobalId( support_id, d_entity_dim );
    const void *tag_data = 0;
    DTK_CHECK_ERROR_CODE( d_moab_mesh->get_moab()->tag_get_by_ptr(
        d_tag, &entity, 1, &tag_data ) );
    const_cast<Scalar *>( static_cast<const Scalar *>( tag_data ) )[dimension] =
        data;
}

//---------------------------------------------------------------------------//
// Finalize a field after writing into it.
template <class Scalar>
void MoabTagField<Scalar>::finalizeAfterWrite()
{
    // Get shared ents.
    moab::Range shared_entities;
    DTK_CHECK_ERROR_CODE( d_moab_mesh->get_shared_entities( -1, shared_entities,
                                                            -1, false, true ) );

    // Exchange the tag.
    DTK_CHECK_ERROR_CODE(
        d_moab_mesh->exchange_tags( d_tag, shared_entities ) );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MOABTAGFIELD_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_MoabTagField_impl.hpp
//---------------------------------------------------------------------------//
