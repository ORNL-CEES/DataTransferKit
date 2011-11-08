//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Transfer_Map.i.hh
 * \author Stuart Slattery
 * \date   Thu Oct 06 15:53:09 2011
 * \brief  Member definitions of class Transfer_Map.
 */
//---------------------------------------------------------------------------//
// $Id: template.i.hh,v 1.4 2008/01/04 22:50:12 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_Transfer_Map_i_hh
#define coupler_Transfer_Map_i_hh

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
void Transfer_Map::add_domain_pair(OrdinateType target_rank, 
				   HandleType source_handle)
{
    target_set.insert(target_rank);
    source_map.insert( Map_Element(target_rank, source_handle) );
}

//---------------------------------------------------------------------------//
// Add a pair to the target map.
void Transfer_Map::add_range_pair(OrdinateType source_rank, 
				  HandleType target_handle)
{
    source_set.insert(source_rank);
    target_map.insert( Map_Element(source_rank, target_handle) );
}

//---------------------------------------------------------------------------//
// Get the number of source handles with a specific target rank.
int Transfer_Map::domain_size(OrdinateType target_rank)
{
    return source_map.count(target_rank);
}
    
//---------------------------------------------------------------------------//
// Get the number of target handles with a specific source rank.
int Transfer_Map::range_size(OrdinateType source_rank)
{
    return target_map.count(source_rank);
}

//---------------------------------------------------------------------------//
// Get the iterator pair for the source domain of a target rank.
Transfer_Map::Map_Pair 
Transfer_Map::source_domain(OrdinateType target_rank)
{
    return source_map.equal_range(target_rank);
}

//---------------------------------------------------------------------------//
// Get the iterator pair for the target range of a source rank.
Transfer_Map::Map_Pair 
Transfer_Map::target_range(OrdinateType source_rank)
{
    return target_map.equal_range(source_rank);
}

//---------------------------------------------------------------------------//
// Return a const_iterator pair to the beginning and end of the source rank
// set. 
Transfer_Map::Set_Pair Transfer_Map::source_set()
{
    
    return Set_Pair( source_set.begin(), source_set.end() );
}

//---------------------------------------------------------------------------//
// Return a const_iterator pair to the beginning and end of the target rank
// set. 
Transfer_Map::Set_Pair Transfer_Map::target_set()
{
    return Set_Pair( target_set.begin() );
}
    
//---------------------------------------------------------------------------//

} // end namespace coupler

#endif // coupler_Transfer_Map_i_hh

//---------------------------------------------------------------------------//
//              end of coupler/Transfer_Map.i.hh
//---------------------------------------------------------------------------//
