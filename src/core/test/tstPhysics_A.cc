//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/test/tstPhysics_A.cc
 * \author stuart
 * \date   Fri Nov 11 09:18:33 2011
 * \brief  
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
#include "Physics_A.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

using physics_A::Physics_A;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void physics_a_test(Parallel_Unit_Test &ut)
{
    // Create a Physics_A instance.
    Physics_A a(MPI_COMM_WORLD,
		0.0, 1.0,
		0.0, 1.0,
		2, 2);

    // Check the communicator.
    UNIT_TEST( a.comm() == MPI_COMM_WORLD );

    // Check the mesh partitioning.
    std::vector<double> x_edges;
    a.get_x_edges(x_edges);

    std::vector<double> y_edges;
    a.get_y_edges(y_edges);

    if (nodes == 1)
    {
	UNIT_TEST( x_edges.size() == 3 );
	UNIT_TEST( x_edges[0] == 0.0 );
	UNIT_TEST( x_edges[1] == 0.5 );
	UNIT_TEST( x_edges[2] == 1.0 );

	UNIT_TEST( x_edges.size() == 3 );
	UNIT_TEST( y_edges[0] == 0.0 );
	UNIT_TEST( y_edges[1] == 0.5 );
	UNIT_TEST( y_edges[2] == 1.0 );
    }
    else if (nodes == 2)
    {
	if (node == 0)
	{
	    UNIT_TEST( x_edges.size() == 2 );
	    UNIT_TEST( x_edges[0] == 0.0 );
	    UNIT_TEST( x_edges[1] == 0.5 );

	    UNIT_TEST( y_edges.size() == 3 );
	    UNIT_TEST( y_edges[0] == 0.0 );
	    UNIT_TEST( y_edges[1] == 0.5 );
	    UNIT_TEST( y_edges[2] == 1.0 );
	}
	else if (node == 1)
	{
	    UNIT_TEST( x_edges.size() == 2 );
	    UNIT_TEST( x_edges[0] == 0.5 );
	    UNIT_TEST( x_edges[1] == 1.0 );

	    UNIT_TEST( y_edges.size() == 3 );
	    UNIT_TEST( y_edges[0] == 0.0 );
	    UNIT_TEST( y_edges[1] == 0.5 );
	    UNIT_TEST( y_edges[2] == 1.0 );
	}
    }
    else if (nodes == 4)
    {
	if (node == 0)
	{
	    UNIT_TEST( x_edges.size() == 2 );
	    UNIT_TEST( x_edges[0] == 0.0 );
	    UNIT_TEST( x_edges[1] == 0.5 );

	    UNIT_TEST( y_edges.size() == 2 );
	    UNIT_TEST( y_edges[0] == 0.0 );
	    UNIT_TEST( y_edges[1] == 0.5 );
	}
	else if (node == 1)
	{
	    UNIT_TEST( x_edges.size() == 2 );
	    UNIT_TEST( x_edges[0] == 0.5 );
	    UNIT_TEST( x_edges[1] == 1.0 );

	    UNIT_TEST( y_edges.size() == 2 );
	    UNIT_TEST( y_edges[0] == 0.0 );
	    UNIT_TEST( y_edges[1] == 0.5 );
	}
	else if (node == 2)
	{
	    UNIT_TEST( x_edges.size() == 2 );
	    UNIT_TEST( x_edges[0] == 0.0 );
	    UNIT_TEST( x_edges[1] == 0.5 );

	    UNIT_TEST( y_edges.size() == 2 );
	    UNIT_TEST( y_edges[0] == 0.5 );
	    UNIT_TEST( y_edges[1] == 1.0 );
	}
	else if (node == 3)
	{
	    UNIT_TEST( x_edges.size() == 2 );
	    UNIT_TEST( x_edges[0] == 0.5 );
	    UNIT_TEST( x_edges[1] == 1.0 );

	    UNIT_TEST( y_edges.size() == 2 );
	    UNIT_TEST( y_edges[0] == 0.5 );
	    UNIT_TEST( y_edges[1] == 1.0 );
	}
    }


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
    
	physics_a_test(ut);
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
        std::cout << "ERROR: While testing tstPhysics_A, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstPhysics_A, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstPhysics_A.cc
//---------------------------------------------------------------------------//
