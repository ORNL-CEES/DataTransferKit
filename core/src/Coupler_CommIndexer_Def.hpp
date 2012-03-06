//---------------------------------------------------------------------------//
/*!
 * \file Coupler_CommIndexer_Def.hpp
 * \author Stuart Slattery
 * \brief CommIndexer definition.
 */
//---------------------------------------------------------------------------//

#ifndef COUPLER_COMMINDEXER_DEF_HPP
#define COUPLER_COMMINDEXER_DEF_HPP

#include <vector>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ENull.hpp>

namespace Coupler
{

/*!
 * \brief Constructor.
 */
template<class T, class Ordinal>
CommIndexer<T,Ordinal>::CommIndexer( RCP_Communicator global_comm, 
				     RCP_Communicator local_comm,
				     RCP_T t_ptr )
{
    int local_rank = local_comm->getRank();
    int global_rank = global_comm->getRank();
    int global_size = global_comm->getSize();

    std::vector<int> t_ids( global_size, 0 );
    if ( t_ptr != Teuchos::null )
    {
    	t_ids[ global_rank ] = 1;
    }
    Teuchos::reduceAll<int,int>( *global_comm,
    				 Teuchos::REDUCE_SUM, 
    				 (int) t_ids.size(),
    				 &t_ids[0], 
    				 &t_ids[0]);

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
   
    std::vector<int>::const_iterator t_it = t_ids.begin();
    typename std::vector<Ordinal>::const_iterator local_it = local_ids.begin();
    typename std::vector<Ordinal>::const_iterator global_it;
    for ( global_it = global_ids.begin(); global_it != global_ids.end();
    	  ++t_it, ++local_it, ++global_it )
    {
    	if ( *t_it )
    	{
    	    d_l2gmap[*local_it] = *global_it;
    	}
    }
}

/*!
 * \brief Destructor.
 */
template<class T, class Ordinal>
CommIndexer<T,Ordinal>::~CommIndexer()
{ /* ... */ }

/*!
 * \brief Given a process id in the local communicator, return the distributed
 * object's process id in the global communicator.
 */
template<class T, class Ordinal>
const Ordinal CommIndexer<T,Ordinal>::l2g( const Ordinal local_id ) const
{
    Ordinal global_id = -1;
    typename IndexMap::const_iterator l2g_pair = d_l2gmap.find( local_id );
    if ( l2g_pair != d_l2gmap.end() )
    {
	global_id = l2g_pair->second;
    }
    return global_id;
}

} // end namespace Coupler

#endif // end COUPLER_COMMINDEXER_DEF_HPP

//---------------------------------------------------------------------------//
// end Coupler_CommIndexer_Def.hpp
//---------------------------------------------------------------------------//
