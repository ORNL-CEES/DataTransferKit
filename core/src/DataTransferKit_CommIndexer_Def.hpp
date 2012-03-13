//---------------------------------------------------------------------------//
/*!
 * \file DataTransferKit_CommIndexer_Def.hpp
 * \author Stuart Slattery
 * \brief CommIndexer definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_COMMINDEXER_DEF_HPP
#define DTK_COMMINDEXER_DEF_HPP

#include <vector>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ENull.hpp>

namespace DataTransferKit
{

/*!
 * \brief Constructor.
 */
template<class Ordinal>
CommIndexer<Ordinal>::CommIndexer( RCP_Communicator global_comm, 
				   RCP_Communicator local_comm )
{
    int local_rank = -1;
    int global_rank = global_comm->getRank();
    int global_size = global_comm->getSize();

    std::vector<int> in_local( global_size, 0 );

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

    std::vector<Ordinal> local_ids( global_size, 0 );
    local_ids[ global_rank ] = local_rank;
    Teuchos::reduceAll<Ordinal,Ordinal>( *global_comm,
					 Teuchos::REDUCE_SUM, 
					 (Ordinal) local_ids.size(),
					 &local_ids[0], 
					 &local_ids[0]);

    std::vector<Ordinal> global_ids( global_size, 0 );
    global_ids[ global_rank ] = global_rank;
    Teuchos::reduceAll<Ordinal,Ordinal>( *global_comm,
					 Teuchos::REDUCE_SUM, 
					 (Ordinal) global_ids.size(),
					 &global_ids[0], 
					 &global_ids[0]);

    std::vector<int>::const_iterator in_local_it = in_local.begin();
    typename std::vector<Ordinal>::const_iterator local_it = 
	local_ids.begin();
    typename std::vector<Ordinal>::const_iterator global_it;
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

/*!
 * \brief Destructor.
 */
template<class Ordinal>
CommIndexer<Ordinal>::~CommIndexer()
{ /* ... */ }

/*!
 * \brief Given a process id in the local communicator, return the distributed
 * object's process id in the global communicator. Return -1 if this local id
 * does not exist in the map.
 */
template<class Ordinal>
const Ordinal CommIndexer<Ordinal>::l2g( const Ordinal local_id ) const
{
    Ordinal global_id = -1;
    typename IndexMap::const_iterator l2g_pair = d_l2gmap.find( local_id );
    if ( l2g_pair != d_l2gmap.end() )
    {
	global_id = l2g_pair->second;
    }
    return global_id;
}

} // end namespace DataTransferKit

#endif // end DTK_COMMINDEXER_DEF_HPP

//---------------------------------------------------------------------------//
// end DataTransferKit_CommIndexer_Def.hpp
//---------------------------------------------------------------------------//
