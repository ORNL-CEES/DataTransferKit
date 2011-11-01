//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstField_DB.cc
 * \author Stuart Slattery
 * \date   Mon Oct 31 12:15:22 2011
 * \brief  Unit tests for the Field_DB class.
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
#include "../SP.hh"
#include "../Field.hh"
#include "../Field_DB.hh"

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

// Database test.
void db_test(Parallel_Unit_Test &ut)
{
    // Build a field database.
    Field_DB<double> field_db;

    // Add the field to the database.
    field_db.add_field("TEST_FIELD");

    // Get the field we made.
    denovo::SP<Field<double> > field = field_db.get_field("TEST_FIELD");

    // Check that the field has its maps.
    UNIT_TEST( field->map_a2b() );
    UNIT_TEST( field->map_b2a() );
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
        
        // add up global passes and fails
        nemesis::global_sum(gpass);
        nemesis::global_sum(gfail);
        ut.numPasses = gpass;
        ut.numFails  = gfail;
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstField_DB, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstField_DB, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstField_DB.cc
//---------------------------------------------------------------------------//
