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
 * \brief DTK_BasicEntityPredicates.cpp
 * \author Stuart R. Slattery
 * \brief Basic entity predicates.
 */
//---------------------------------------------------------------------------//

#include "DTK_BasicEntityPredicates.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Block predicate.
bool BlockPredicate::operator()( Entity entity )
{
    bool found_in_block = false;
    Teuchos::Array<int>::const_iterator block_it;
    for ( block_it = d_block_ids.begin(); block_it != d_block_ids.end();
          ++block_it )
    {
        if ( entity.inBlock( *block_it ) )
        {
            found_in_block = true;
            break;
        }
    }
    return found_in_block;
}

//---------------------------------------------------------------------------//
// Boundary predicate.
bool BoundaryPredicate::operator()( Entity entity )
{
    bool found_on_boundary = false;
    Teuchos::Array<int>::const_iterator boundary_it;
    for ( boundary_it = d_boundary_ids.begin();
          boundary_it != d_boundary_ids.end(); ++boundary_it )
    {
        if ( entity.onBoundary( *boundary_it ) )
        {
            found_on_boundary = true;
            break;
        }
    }
    return found_on_boundary;
}

//---------------------------------------------------------------------------//
// LocalEntity predicate.
bool LocalEntityPredicate::operator()( Entity entity )
{
    return ( d_my_rank == entity.ownerRank() );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_BasicEntityPredicates.cpp
//---------------------------------------------------------------------------//
