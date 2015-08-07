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
 * \brief DTK_MoabEntityPredicates.cpp
 * \author Stuart R. Slattery
 * \brief Moab entity predicates.
 */
//---------------------------------------------------------------------------//

#include "DTK_MoabEntityPredicates.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Single mesh set constructor.
MoabMeshSetPredicate::MoabMeshSetPredicate( 
    const moab::EntityHandle& mesh_set,
    const Teuchos::RCP<MoabMeshSetIndexer>& set_indexer )
    : d_set_ids( 1, set_indexer->getIndexFromMeshSet(mesh_set) )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Vector Constructor.
MoabMeshSetPredicate::MoabMeshSetPredicate( 
    const std::vector<moab::EntityHandle>& mesh_sets,
    const Teuchos::RCP<MoabMeshSetIndexer>& set_indexer )
    : d_set_ids( mesh_sets.size(), -1 )
{
    int num_sets = d_set_ids.size();
    for ( int i = 0; i < num_sets; ++i )
    {
	d_set_ids[i] = set_indexer->getIndexFromMeshSet( mesh_sets[i] );
    }
}

//---------------------------------------------------------------------------//
// Functor.
bool MoabMeshSetPredicate::operator()( Entity entity ) 
{ 
    for ( auto set_it = d_set_ids.begin();
	  set_it != d_set_ids.end();
	  ++set_it )
    {
	if ( !entity.inBlock(*set_it) )
	{
	    return false;
	}
    }
    return true;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MoabEntityPredicates.cpp
//---------------------------------------------------------------------------//
