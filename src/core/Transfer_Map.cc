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

namespace dtransfer
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
// Add a rank/index to the map.
void Transfer_Map::add_pair(int rank, int index)
{
    element_rank.push_back(rank);
    element_index.push_back(index);
    Check( index.size() == rank.size() );
}

//---------------------------------------------------------------------------//

} // end namespace dtransfer

//---------------------------------------------------------------------------//
//                 end of Transfer_Map.cc
//---------------------------------------------------------------------------//
