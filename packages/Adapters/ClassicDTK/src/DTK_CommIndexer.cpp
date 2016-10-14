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
 * \file DTK_CommIndexer.cpp
 * \author Stuart Slattery
 * \brief CommIndexer definition.
 */
//---------------------------------------------------------------------------//

#include "DTK_CommIndexer.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_as.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
CommIndexer::CommIndexer() { /* ... */}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor. Given a global communicator and a local communicator,
 * index the local communicator process ranks into the global communicator
 * process ranks.
 *
 * \param global_comm The global communicator.
 *
 * \param local_comm The local communicator.
 */
CommIndexer::CommIndexer( Teuchos::RCP<const Teuchos::Comm<int>> global_comm,
                          Teuchos::RCP<const Teuchos::Comm<int>> local_comm )
{
    // Set whether or not the indexer will be valid on this rank.
    d_is_valid = Teuchos::nonnull( local_comm );

    // Get my rank in the local communicator.
    int local_rank = d_is_valid ? local_comm->getRank() : -1;

    // Gather everyone's rank in the local communicator.
    int global_size = global_comm->getSize();
    Teuchos::Array<int> local_ids( global_size, 0 );
    Teuchos::gatherAll<int, int>( *global_comm, 1, &local_rank,
                                  local_ids.size(), local_ids.getRawPtr() );

    // Map the local communicator to the global communicator.
    for ( int i = 0; i < global_size; ++i )
    {
        if ( local_ids[i] >= 0 )
        {
            d_l2gmap[local_ids[i]] = i;
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
CommIndexer::~CommIndexer() { /* ... */}

//---------------------------------------------------------------------------//
/*!
 * \brief Given a process id in the local communicator, return the distributed
 * object's process id in the global communicator.
 *
 * \param local_id The local communicator process rank.
 *
 * \return The global process rank. Return -1 if this local id does not exist
 * in the map.
 */
int CommIndexer::l2g( const int local_id ) const
{
    std::unordered_map<int, int>::const_iterator l2g_pair =
        d_l2gmap.find( local_id );
    return ( l2g_pair != d_l2gmap.end() ) ? l2g_pair->second : -1;
}

//---------------------------------------------------------------------------//
//! Return the size of the local to global map.
int CommIndexer::size() const { return d_l2gmap.size(); }

//---------------------------------------------------------------------------//
// Return true if the indexer is valid on this process (local_comm is
// nonnull).
bool CommIndexer::isValid() const { return d_is_valid; }

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_CommIndexer.cpp
//---------------------------------------------------------------------------//
