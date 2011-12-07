//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstTransfer_Map.cc
 * \author Stuart Slattery
 * \date   Mon Oct 31 12:15:31 2011
 * \brief  Unit test for the transfer map.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <map>
#include <set>
#include <iterator>

#include "Teuchos_UnitTestHarness.hpp"

#include "../Transfer_Map.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

namespace coupler {

// Test the Transfer_Map data structure.
TEUCHOS_UNIT_TEST( Transfer_Map, transfermapStorage )
{
    // Useful typedefs.
    typedef Transfer_Map::Map_Iterator    Map_Iterator;
    typedef Transfer_Map::Map_Pair        Map_Pair;
    typedef Transfer_Map::Set_Iterator    Set_Iterator;
    typedef Transfer_Map::Set_Pair        Set_Pair;

    // Make a map object with integer handles and ranks.
    Transfer_Map tmap;

    // Add domain pairs.
    tmap.add_domain_pair(1, 21);
    tmap.add_domain_pair(1, 32);
    tmap.add_domain_pair(1, 12);
    tmap.add_domain_pair(3, 2);
    tmap.add_domain_pair(4, 22);
    tmap.add_domain_pair(2, 54);
    tmap.add_domain_pair(3, 5);

    // Check the size of the domains.
    TEST_ASSERT( tmap.domain_size(1) == 3 );
    TEST_ASSERT( tmap.domain_size(2) == 1 );
    TEST_ASSERT( tmap.domain_size(3) == 2 );
    TEST_ASSERT( tmap.domain_size(4) == 1 );

    // Check the contents of the domains.
    Map_Pair domain_1 = tmap.domain(1);
    TEST_ASSERT( std::distance(domain_1.first, domain_1.second) == 3 );
    Map_Iterator domain_it_1 = domain_1.first;
    TEST_ASSERT( domain_it_1->second == 21 );
    ++domain_it_1;
    TEST_ASSERT( domain_it_1->second == 32 );
    ++domain_it_1;
    TEST_ASSERT( domain_it_1->second == 12 );
    
    Map_Pair domain_2 = tmap.domain(2);
    TEST_ASSERT( std::distance(domain_2.first, domain_2.second) == 1 );
    Map_Iterator domain_it_2 = domain_2.first;
    TEST_ASSERT( domain_it_2->second == 54 );

    Map_Pair domain_3 = tmap.domain(3);
    TEST_ASSERT( std::distance(domain_3.first, domain_3.second) == 2 );
    Map_Iterator domain_it_3 = domain_3.first;
    TEST_ASSERT( domain_it_3->second == 2 );
    ++domain_it_3;
    TEST_ASSERT( domain_it_3->second == 5 );

    Map_Pair domain_4 = tmap.domain(4);
    TEST_ASSERT( std::distance(domain_4.first, domain_4.second) == 1 );
    Map_Iterator domain_it_4 = domain_4.first;
    TEST_ASSERT( domain_it_4->second == 22 );

    // Add range pairs.
    tmap.add_range_pair(8, 9);
    tmap.add_range_pair(32, 1);
    tmap.add_range_pair(8, 4);
    tmap.add_range_pair(8, 34);
    tmap.add_range_pair(32, 3);

    // Check the size of the ranges.
    TEST_ASSERT( tmap.range_size(8) == 3 );
    TEST_ASSERT( tmap.range_size(32) == 2 );

    // Check the contents of the ranges.
    Map_Pair range_8 = tmap.range(8);
    TEST_ASSERT( std::distance(range_8.first,range_8.second) == 3 );
    Map_Iterator range_it_8 = range_8.first;
    TEST_ASSERT( range_it_8->second == 9 );
    ++range_it_8;
    TEST_ASSERT( range_it_8->second == 4 );
    ++range_it_8;
    TEST_ASSERT( range_it_8->second == 34 );

    Map_Pair range_32 = tmap.range(32);
    TEST_ASSERT( std::distance(range_32.first,range_32.second) == 2 );
    Map_Iterator range_it_32 = range_32.first;
    TEST_ASSERT( range_it_32->second == 1 );
    ++range_it_32;
    TEST_ASSERT( range_it_32->second == 3 );

    // Check the source and target sets.
    Set_Pair sources = tmap.sources();
    TEST_ASSERT( std::distance(sources.first, sources.second) == 2 );
    Set_Iterator source_it = sources.first;
    TEST_ASSERT( *source_it == 8 );
    ++source_it;
    TEST_ASSERT( *source_it == 32 );

    Set_Pair targets = tmap.targets();
    TEST_ASSERT( std::distance(targets.first, targets.second) == 4 );
    Set_Iterator target_it = targets.first;
    TEST_ASSERT( *target_it == 1 );
    ++target_it;
    TEST_ASSERT( *target_it == 2 );
    ++target_it;
    TEST_ASSERT( *target_it == 3 );
    ++target_it;
    TEST_ASSERT( *target_it == 4 );
}

} // end namespace coupler

//---------------------------------------------------------------------------//
//                        end of tstTransfer_Map.cc
//---------------------------------------------------------------------------//
