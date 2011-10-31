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

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "comm/global.hh"
#include "comm/Parallel_Unit_Test.hh"
#include "release/Release.hh"
#include "../Transfer_Map.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

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
    typedef typename std::multimap<int,int>::const_iterator   Map_Iterator;
    typedef std::pair<Map_Iterator, Map_Iterator>             Iterator_Pair;
    typedef typename std::set<int>::const_iterator            Set_Iterator;

    // Make a map object with integer handles and ranks.
    Transfer_Map<int,int> tmap;

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
    Iterator_Pair domain_1 = tmap.domain(1);
    UNIT_TEST( std::distance(domain_1.first(), domain_1.second()) == 3 );
    UNIT_TEST( *(domain_1.first().second()) == 12 );
    UNIT_TEST( *(domain_1.first().second()+1) == 21 );
    UNIT_TEST( *(domain_1.first().second()+2) == 32 );

    Iterator_Pair domain_2 = tmap.domain(2);
    UNIT_TEST( std::distance(domain_2.first().second(), domain_2.second()) == 1 );
    UNIT_TEST( *(domain_2.first().second()) == 54 );

    Iterator_Pair domain_3 = tmap.domain(3);
    UNIT_TEST( std::distance(domain_3.first().second(), domain_3.second()) == 2 );
    UNIT_TEST( *(domain_3.first().second()) == 2 );
    UNIT_TEST( *(domain_3.first().second()+1) == 5 );

    Iterator_Pair domain_4 = tmap.domain(4);
    UNIT_TEST( std::distance(domain_4.first().second(), domain_4.second()) == 1 );
    UNIT_TEST( *(domain_4.first().second()) == 22 );

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
    Iterator_Pair range_8 = tmap.range(8);
    UNIT_TEST( std::distance(range_8.first().second(),range_8.second()) == 3 );
    UNIT_TEST( *(range_8.first().second()) == 4 );
    UNIT_TEST( *(range_8.first().second()+1) == 9 );
    UNIT_TEST( *(range_8.first().second()+2) == 34 );

    Iterator_Pair range_32 = tmap.range(32);
    UNIT_TEST( std::distance(range_32.first().second(),range_32.second()) == 2 );
    UNIT_TEST( *(range_32.first().second()) == 1 );
    UNIT_TEST( *(range_32.first().second()+1) == 3 );

    // Check the source and target sets.
    Set_Iterator source_begin = tmap.source_set_begin();
    Set_Iterator source_end = tmap.source_set_end();
    UNIT_TEST( std::distance(source_begin, source_end) == 2 );
    UNIT_TEST( *source_begin == 8 );
    UNIT_TEST( *(source_begin+1) == 32 );

    Set_Iterator target_begin = tmap.target_set_begin();
    Set_Iterator target_end = tmap.target_set_end();
    UNIT_TEST( std::distance(target_begin, target_end) == 4 );
    UNIT_TEST( *target_begin == 1 );
    UNIT_TEST( *(target_begin+1) == 2 );
    UNIT_TEST( *(target_begin+2) == 3 );
    UNIT_TEST( *(target_begin+3) == 4 );

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
    Parallel_Unit_Test ut(argc, argv, release);

    node  = nemesis::node();
    nodes = nemesis::nodes();
    
    try
    {
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;

	// The map only needs to be tested in serial
	if (node == 0)
	{
	    map_test(ut);
	}

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
