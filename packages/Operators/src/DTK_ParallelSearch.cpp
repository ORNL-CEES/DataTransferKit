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
#include "DTK_EntitySequence.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
ParallelSearch::ParallelSearch( 
    const Teuchos::RCP<FunctionSpace>& domain_space,
    const std::function<bool(Entity)>& domain_predicate,
    const EntityType domain_entity_type,
    const Teuchos::RCP<FunctionSpace>& range_space,
    const std::function<bool(Entity)>& range_predicate
    const EntityType range_entity_type )
{
    // COARSE GLOBAL SEARCH
    // Determine which parallel regions for the domain intersect which
    // parallel regions of the range.

    // Determine which centroids in the local range parallel region
    // will be sent to which domain parallel regions.

    // Send the range centroids.

    // COARSE LOCAL SEARCH
    // Construct a search tree from the local domain centroids.

    // Find the set of nearest domain entities to each range centroid.

    // FINE LOCAL SEARCH
    // Map the range centroid to the reference frame of the nearest domain
    // entities to determine point inclusion.

    // Store the parametric coordinates.

    // Add the entry to the graph.
}

//---------------------------------------------------------------------------//
// Destructor.
ParallelSearch::~ParallelSearch()
{ /* ... */ }

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_PARALLELSEARCH_HPP

//---------------------------------------------------------------------------//
// end DTK_ParallelSearch.hpp
//---------------------------------------------------------------------------//
