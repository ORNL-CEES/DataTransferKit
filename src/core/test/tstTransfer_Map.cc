//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstTransfer_Map.cc
 * \author Stuart Slattery
 * \date   Mon Oct 31 12:15:31 2011
 * \brief  Unit test for the transfer map.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template_c4_test.cc,v 1.7 2008/01/02 22:50:26 9te Exp $
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <map>
#include <set>
#include <iterator>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "comm/global.hh"
#include "comm/Parallel_Unit_Test.hh"
#include "release/Release.hh"
#include "../Transfer_Map.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

using coupler::Transfer_Map;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void map_test(Parallel_Unit_Test &ut)
{
    // Useful typedefs.
    typedef typename Transfer_Map::Map_Iterator    Map_Iterator;
    typedef typename Transfer_Map::Map_Pair        Map_Pair;
    typedef typename Transfer_Map::Set_Iterator    Set_Iterator;
    typedef typename Transfer_Map::Set_Pair        Set_Pair;

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
    UNIT_TEST( tmap.domain_size(1) == 3 );
    UNIT_TEST( tmap.domain_size(2) == 1 );
    UNIT_TEST( tmap.domain_size(3) == 2 );
    UNIT_TEST( tmap.domain_size(4) == 1 );

    // Check the contents of the domains.
    Map_Pair domain_1 = tmap.domain(1);
    UNIT_TEST( std::distance(domain_1.first, domain_1.second) == 3 );
    Map_Iterator domain_it_1 = domain_1.first;
    UNIT_TEST( domain_it_1->second == 21 );
    ++domain_it_1;
    UNIT_TEST( domain_it_1->second == 32 );
    ++domain_it_1;
    UNIT_TEST( domain_it_1->second == 12 );
    
    Map_Pair domain_2 = tmap.domain(2);
    UNIT_TEST( std::distance(domain_2.first, domain_2.second) == 1 );
    Map_Iterator domain_it_2 = domain_2.first;
    UNIT_TEST( domain_it_2->second == 54 );

    Map_Pair domain_3 = tmap.domain(3);
    UNIT_TEST( std::distance(domain_3.first, domain_3.second) == 2 );
    Map_Iterator domain_it_3 = domain_3.first;
    UNIT_TEST( domain_it_3->second == 2 );
    ++domain_it_3;
    UNIT_TEST( domain_it_3->second == 5 );

    Map_Pair domain_4 = tmap.domain(4);
    UNIT_TEST( std::distance(domain_4.first, domain_4.second) == 1 );
    Map_Iterator domain_it_4 = domain_4.first;
    UNIT_TEST( domain_it_4->second == 22 );

    // Add range pairs.
    tmap.add_range_pair(8, 9);
    tmap.add_range_pair(32, 1);
    tmap.add_range_pair(8, 4);
    tmap.add_range_pair(8, 34);
    tmap.add_range_pair(32, 3);

    // Check the size of the ranges.
    UNIT_TEST( tmap.range_size(8) == 3 );
    UNIT_TEST( tmap.range_size(32) == 2 );

    // Check the contents of the ranges.
    Map_Pair range_8 = tmap.range(8);
    UNIT_TEST( std::distance(range_8.first,range_8.second) == 3 );
    Map_Iterator range_it_8 = range_8.first;
    UNIT_TEST( range_it_8->second == 9 );
    ++range_it_8;
    UNIT_TEST( range_it_8->second == 4 );
    ++range_it_8;
    UNIT_TEST( range_it_8->second == 34 );

    Map_Pair range_32 = tmap.range(32);
    UNIT_TEST( std::distance(range_32.first,range_32.second) == 2 );
    Map_Iterator range_it_32 = range_32.first;
    UNIT_TEST( range_it_32->second == 1 );
    ++range_it_32;
    UNIT_TEST( range_it_32->second == 3 );

    // Check the source and target sets.
    Set_Pair sources = tmap.sources();
    UNIT_TEST( std::distance(sources.first, sources.second) == 2 );
    Set_Iterator source_it = sources.first;
    UNIT_TEST( *source_it == 8 );
    ++source_it;
    UNIT_TEST( *source_it == 32 );

    Set_Pair targets = tmap.targets();
    UNIT_TEST( std::distance(targets.first, targets.second) == 4 );
    Set_Iterator target_it = targets.first;
    UNIT_TEST( *target_it == 1 );
    ++target_it;
    UNIT_TEST( *target_it == 2 );
    ++target_it;
    UNIT_TEST( *target_it == 3 );
    ++target_it;
    UNIT_TEST( *target_it == 4 );

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "Transfer_Map test passes on " << node;
        ut.passes( m.str() );
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Parallel_Unit_Test ut(argc, argv, coupler::release);

    node  = nemesis::node();
    nodes = nemesis::nodes();
    
    try
    {
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;

	map_test(ut);
	gpass += ut.numPasses;
	gfail += ut.numFails;
	ut.reset();
    
        // add up global passes and fails
        nemesis::global_sum(gpass);
        nemesis::global_sum(gfail);
	ut.numPasses = gpass;
        ut.numFails  = gfail;
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstTransfer_Map, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstTransfer_Map, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstTransfer_Map.cc
//---------------------------------------------------------------------------//
