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

#ifndef core_Transfer_Map_i_hh
#define core_Transfer_Map_i_hh

#include "harness/DBC.hh"

namespace coupler
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Transfer_Map::Transfer_Map()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Transfer_Map::~Transfer_Map()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Add a pair to the source map.The source handle corresponding to an
 * entity in the local domain correlates to the range owned by the target
 * rank. 
 * \param target_rank The MPI ordinate of the target.
 * \param handle The local source handle associated with the target_rank.
 */
void Transfer_Map::add_domain_pair(OrdinateType target_rank, 
				   HandleType handle)
{
    target_set.insert(target_rank);
    source_map.insert( Map_Element(target_rank, source_handle) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add a pair to the target map. The target handle corresponding to an
 * entity in the local range correlates to the domain owned by the source
 * rank. 
 * \param source_rank The MPI ordinate of the source.
 * \param handle The local target handle associated with the source_rank.
 */
void Transfer_Map::add_range_pair(OrdinateType source_rank, 
				  HandleType handle)
{
    source_set.insert(source_rank);
    target_map.insert( Map_Element(source_rank, target_handle) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the number of source handles with a specific target rank. This
 * is the size of the local domain correlating to the target rank. 
 * \param target_rank The target rank for which the local domain size is
 * desired. 
 * \return The size of the local domain owned by the target rank.
 */
int Transfer_Map::domain_size(OrdinateType target_rank)
{
    return source_map.count(target_rank);
}
    
//---------------------------------------------------------------------------//
/*!
 * \brief Get the number of target handles with a specific source rank. This
 * is the size of the local range correlating to the source rank. 
 * \param source_rank The source rank for which the local range size is
 * desired. 
 * \return The size of the local range owned by the source rank.
 */
int Transfer_Map::range_size(OrdinateType source_rank)
{
    return target_map.count(source_rank);
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Get the iterator pair for the source domain of a target rank.This
 * correlates to the local source handles that exist in the range of the
 * target rank. 
 * \param target_rank The target rank for which the local source domain is
 * desired. 
 * \return An iterator pair over the local domain for the given target rank.
 */
Transfer_Map::Map_Pair Transfer_Map::source_domain(OrdinateType target_rank)
{
    return source_map.equal_range(target_rank);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the iterator pair for the target range of a source rank. This
 * correlates to the local target handles that exist in the domain of the
 * source rank. 
 * \param source_rank The source rank for which the local target range is
 * desired. 
 * \return An iterator pair over the local range for the give source rank.
 */
Transfer_Map::Map_Pair Transfer_Map::target_range(OrdinateType source_rank)
{
    return target_map.equal_range(source_rank);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a const_iterator pair to the beginning and end of the source
 * rank set.
 */
Transfer_Map::Set_Pair Transfer_Map::source_set()
{
    
    return Set_Pair( source_set.begin(), source_set.end() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a const_iterator pair to the beginning and end of the target
 * rank set.
 */
Transfer_Map::Set_Pair Transfer_Map::target_set()
{
    return Set_Pair( target_set.begin() );
}
    
//---------------------------------------------------------------------------//

} // end namespace coupler

#endif // core_Transfer_Map_i_hh

//---------------------------------------------------------------------------//
//              end of core/Transfer_Map.i.hh
//---------------------------------------------------------------------------//
