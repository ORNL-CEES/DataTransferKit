//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   src/Transfer_Map.cc
 * \author Stuart Slattery
 * \date   Thu Oct 06 15:53:09 2011
 * \brief  Transfer_Map member definitions.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "Transfer_Map.hh"
#include "harness/DBC.hh"

namespace coupler
{

//---------------------------------------------------------------------------//
// Constructor.
Transfer_Map::Transfer_Map()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
Transfer_Map::~Transfer_Map()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Add a pair to the source map.
void Transfer_Map::add_domain_pair(int target_rank, int source_handle)
{
    target_set.insert(target_rank);
    source_map.insert( Map_Pair(target_rank, source_handle) );
}

//---------------------------------------------------------------------------//
// Add a pair to the target map.
void Transfer_Map::add_range_pair(int source_rank, int target_handle)
{
    source_set.insert(source_rank);
    target_map.insert( Map_Pair(source_rank, target_handle) );
}

//---------------------------------------------------------------------------//
// Get the number of source handles with a specific target rank.
int Transfer_Map::domain_size(int target_rank)
{
    return source_map.count(target_rank);
}
    
//---------------------------------------------------------------------------//
// Get the number of target handles with a specific source rank.
int Transfer_Map::range_size(int source_rank)
{
    return target_map.count(source_rank);
}


//---------------------------------------------------------------------------//
// Get the iterator pair for the source domain of a target rank.
Transfer_Map::Iterator_Pair Transfer_Map::source_domain(int target_rank)
{
    return source_map.equal_range(target_rank);
}

//---------------------------------------------------------------------------//
// Get the iterator pair for the target range of a source rank.
Transfer_Map::Iterator_Pair Transfer_Map::target_range(int source_rank)
{
    return target_map.equal_range(source_rank);
}

//---------------------------------------------------------------------------//
// Return a const_iterator to the beginning of the source rank set.
Transfer_Map::Set_Iterator Transfer_Map::source_set_begin()
{
    return source_set.begin();
}

//---------------------------------------------------------------------------//
// Return a const_iterator to the end of the source rank set.
Transfer_Map::Set_Iterator Transfer_Map::source_set_end()
{
    return source_set.end();
}

//---------------------------------------------------------------------------//
// Return a const_iterator to the beginning of the target rank set.
Transfer_Map::Set_Iterator Transfer_Map::target_set_begin()
{
    return target_set.begin();
}
    
//---------------------------------------------------------------------------//
// Return a const_iterator to the end of the target rank set.
Transfer_Map::Set_Iterator Transfer_Map::target_set_end()
{
    return target_set.end()
}

//---------------------------------------------------------------------------//

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Transfer_Map.cc
//---------------------------------------------------------------------------//
