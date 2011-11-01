//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstField.cc
 * \author Stuart Slattery
 * \date   Mon Oct 31 12:15:17 2011
 * \brief  Unit tests for the Field class.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template_c4_test.cc,v 1.7 2008/01/02 22:50:26 9te Exp $
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "comm/global.hh"
#include "comm/Parallel_Unit_Test.hh"
#include "release/Release.hh"
#include "../Field.hh"

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

// Test the constructor to check map initialization.
void constructor_test(Parallel_Unit_Test &ut)
{
    // Build a field.
    Field<double> dbl_field;

    // Check that the maps have been constructed.
    UNIT_TEST( dbl_field.map_a2b() );
    UNIT_TEST( dbl_field.map_b2a() );
}

// Test that we can manipulate the map through the field.
void map_test(Parallel_Unit_Test &ut)
{
    // build a field
    Field<double> dbl_field;
    
    // Add a domain pair to the map.
    dbl_field.map_a2b()->add_domain_pair(2, 32);

    // Check the map to see that the domain pair was properly added.
    UNIT_TEST( dbl_field.map_a2b()->domain_size(2) == 1 );
    UNIT_TEST( dbl_field.map_a2b()->domain(2).first().second() == 32 );
    UNIT_TEST( *(dbl_field.map_a2b()->target_set_begin()) == 2 );
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

	constructor_test(ut);
	gpass += ut.numPasses;
        gfail += ut.numFails;
        ut.reset();

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
        std::cout << "ERROR: While testing tstField, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstField, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstField.cc
//---------------------------------------------------------------------------//
