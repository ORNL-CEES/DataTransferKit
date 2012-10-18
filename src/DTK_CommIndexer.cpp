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

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
CommIndexer::CommIndexer()
{ /* ... */ }

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
CommIndexer::CommIndexer( RCP_Comm global_comm, RCP_Comm local_comm )
{
    int local_rank = -1;
    int global_rank = global_comm->getRank();
    int global_size = global_comm->getSize();

    Teuchos::Array<int> in_local( global_size, 0 );
    Teuchos::Array<int> in_local_copy( global_size, 0 );

    if ( !local_comm.is_null() )
    {
	local_rank = local_comm->getRank();
    	in_local[ global_rank ] = 1;
    }
    global_comm->barrier();

    Teuchos::reduceAll<int,int>( *global_comm,
    				 Teuchos::REDUCE_SUM, 
    				 Teuchos::as<int>(in_local.size()),
    				 &in_local[0], 
    				 &in_local_copy[0]);

    Teuchos::Array<int> local_ids( global_size, 0 );
    Teuchos::Array<int> local_ids_copy( global_size, 0 );
    local_ids[ global_rank ] = local_rank;
    Teuchos::reduceAll<int,int>( *global_comm,
				 Teuchos::REDUCE_SUM, 
				 Teuchos::as<int>(local_ids.size()),
				 &local_ids[0], 
				 &local_ids_copy[0]);

    Teuchos::Array<int> global_ids( global_size, 0 );
    Teuchos::Array<int> global_ids_copy( global_size, 0 );
    global_ids[ global_rank ] = global_rank;
    Teuchos::reduceAll<int,int>( *global_comm,
				 Teuchos::REDUCE_SUM, 
				 Teuchos::as<int>(global_ids.size()),
				 &global_ids[0], 
				 &global_ids_copy[0]);

    Teuchos::Array<int>::const_iterator in_local_copy_it = in_local_copy.begin();
    Teuchos::Array<int>::const_iterator local_it = 
	local_ids_copy.begin();
    Teuchos::Array<int>::const_iterator global_it;
    for ( global_it = global_ids_copy.begin(); global_it != global_ids_copy.end();
    	  ++in_local_copy_it, ++local_it, ++global_it )
    {
    	if ( *in_local_copy_it )
    	{
    	    d_l2gmap[*local_it] = *global_it;
    	}
    }

    global_comm->barrier();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
CommIndexer::~CommIndexer()
{ /* ... */ }

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
const int CommIndexer::l2g( const int local_id ) const
{
    int global_id = -1;
    IndexMap::const_iterator l2g_pair = d_l2gmap.find( local_id );
    if ( l2g_pair != d_l2gmap.end() )
    {
	global_id = l2g_pair->second;
    }
    return global_id;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_CommIndexer.cpp
//---------------------------------------------------------------------------//
