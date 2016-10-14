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
 * \brief DTK_STKMeshEntityPredicates.cpp
 * \author Stuart R. Slattery
 * \brief STK mesh entity predicates.
 */
//---------------------------------------------------------------------------//

#include "DTK_STKMeshEntityPredicates.hpp"

#include <stk_mesh/base/MetaData.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Part predicate.
//---------------------------------------------------------------------------//
// Functor.
bool STKPartPredicate::operator()( Entity entity )
{
    for ( auto part_it = b_part_ids.begin(); part_it != b_part_ids.end();
          ++part_it )
    {
        if ( !entity.inBlock( *part_it ) )
        {
            return false;
        }
    }
    return true;
}

//---------------------------------------------------------------------------//
// Part name predicate.
//---------------------------------------------------------------------------//
STKPartNamePredicate::STKPartNamePredicate(
    const Teuchos::Array<std::string> &part_names,
    const Teuchos::RCP<stk::mesh::BulkData> &bulk_data )
{
    stk::mesh::Part *part = 0;
    for ( auto name_it = part_names.begin(); name_it != part_names.end();
          ++name_it )
    {
        part = bulk_data->mesh_meta_data().get_part( *name_it );
        this->b_part_ids.push_back( part->mesh_meta_data_ordinal() );
    }
}

//---------------------------------------------------------------------------//
// Part Vector predicate.
//---------------------------------------------------------------------------//
STKPartVectorPredicate::STKPartVectorPredicate(
    const stk::mesh::PartVector &parts )
{
    for ( auto part_it = parts.begin(); part_it != parts.end(); ++part_it )
    {
        this->b_part_ids.push_back( ( *part_it )->mesh_meta_data_ordinal() );
    }
}

//---------------------------------------------------------------------------//
// Selector predicate.
//---------------------------------------------------------------------------//
STKSelectorPredicate::STKSelectorPredicate(
    const stk::mesh::Selector &selector )
{
    stk::mesh::PartVector parts;
    selector.get_parts( parts );
    for ( auto part_it = parts.begin(); part_it != parts.end(); ++part_it )
    {
        this->b_part_ids.push_back( ( *part_it )->mesh_meta_data_ordinal() );
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_STKMeshEntityPredicates.cpp
//---------------------------------------------------------------------------//
