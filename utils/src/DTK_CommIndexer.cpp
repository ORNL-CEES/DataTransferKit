//---------------------------------------------------------------------------//
/*!
 * \file DTK_CommIndexer.cpp
 * \author Stuart Slattery
 * \brief CommIndexer definition.
 */
//---------------------------------------------------------------------------//

#include "DTK_CommIndexer.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_Array.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
CommIndexer::CommIndexer( RCP_Comm global_comm, RCP_Comm local_comm )
{
    int local_rank = -1;
    int global_rank = global_comm->getRank();
    int global_size = global_comm->getSize();

    Teuchos::Array<int> in_local( global_size, 0 );

    if ( local_comm != Teuchos::null )
    {
	local_rank = local_comm->getRank();
    	in_local[ global_rank ] = 1;
    }
    global_comm->barrier();

    Teuchos::reduceAll<int,int>( *global_comm,
    				 Teuchos::REDUCE_SUM, 
    				 (int) in_local.size(),
    				 &in_local[0], 
    				 &in_local[0]);

    Teuchos::Array<int> local_ids( global_size, 0 );
    local_ids[ global_rank ] = local_rank;
    Teuchos::reduceAll<int,int>( *global_comm,
				 Teuchos::REDUCE_SUM, 
				 (int) local_ids.size(),
				 &local_ids[0], 
				 &local_ids[0]);

    Teuchos::Array<int> global_ids( global_size, 0 );
    global_ids[ global_rank ] = global_rank;
    Teuchos::reduceAll<int,int>( *global_comm,
				 Teuchos::REDUCE_SUM, 
				 (int) global_ids.size(),
				 &global_ids[0], 
				 &global_ids[0]);

    Teuchos::Array<int>::const_iterator in_local_it = in_local.begin();
    typename Teuchos::Array<int>::const_iterator local_it = 
	local_ids.begin();
    typename Teuchos::Array<int>::const_iterator global_it;
    for ( global_it = global_ids.begin(); global_it != global_ids.end();
    	  ++in_local_it, ++local_it, ++global_it )
    {
    	if ( *in_local_it )
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
 * object's process id in the global communicator. Return -1 if this local id
 * does not exist in the map.
 */
const int CommIndexer::l2g( const int local_id ) const
{
    int global_id = -1;
    typename IndexMap::const_iterator l2g_pair = d_l2gmap.find( local_id );
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
