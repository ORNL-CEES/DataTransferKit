//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/test/tstMap_Node.cc
 * \author Stuart R. Slattery
 * \date   Thu May 26 16:02:19 2011
 * \brief  Map_Node class test
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
#include "mesh_type/Point.hh"
#include "utils/SP.hh"
#include "fields/Field_View.hh"
#include "../Map_Node.hh"

using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);


//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
// create an interleaved set of field points
template<typename DataType>
std::vector<DataType> make_field_pts()
{
    std::vector<DataType> pts;
    DataType counter = 0.5;
    for (unsigned int n = 0; n < 10; ++n)
    {
        for (unsigned int i = 0; i < 3; ++i) 
        {
            pts.push_back(counter);
        }

        counter += 1.0;
    }
    return pts;
}

//---------------------------------------------------------------------------//
// make an arbitrary data field
template<typename DataType>
std::vector<DataType> make_data_field()
{
    std::vector<DataType> data;
    DataType counter = 0.5;
    for (unsigned int n = 0; n < 10; ++n) 
    {
        data.push_back(counter);
        counter += 1.0;
    }
    return data;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// Tests point-based map node
void serial_node_point_test(Parallel_Unit_Test &ut)
{ 
    typedef double                                      DataType;
    typedef denovo::Point<DataType>                     Point_t;
    typedef coupler::Map_Node<int, DataType, Point_t>   Map_Node_t;
    typedef std::vector<DataType>                       Vec_DataType;

    if (nodes > 1)
        return;

    // make data field
    Vec_DataType data_field = make_data_field<DataType>();

    // make field points
    Vec_DataType pts = make_field_pts<DataType>();

    // make a vector to hold the nodes
    std::vector<Map_Node_t> nodes;

    // fill the node vector
    unsigned int partition = 0;
    for (unsigned int n = 0; n < 10; ++n)
    {
        // make a new point
        Point_t p(pts[3*n], pts[3*n+1], pts[3*n+2]);

        // make a new node
        Map_Node_t node(n, partition, p);

        // verify that power, temperature and volume do not currently exist
        UNIT_TEST( !node.exists("power") );
        UNIT_TEST( !node.exists("temperature") );
        UNIT_TEST( !node.exists("volume") );

        // assign some data to the node
        node.add("power", data_field[n]);
        node.add("temperature", 2*data_field[n]);
        node.add("volume", 3*data_field[n]);

        // verify that power, temperature and volume do currently exist
        UNIT_TEST( node.exists("power") );
        UNIT_TEST( node.exists("temperature") );
        UNIT_TEST( node.exists("volume") );

        // add it to the node vector
        nodes.push_back(node);

        // update the psuedo partition count
        partition += 2;
    }

    // verify the contents of the nodes
    DataType counter = 0.5;
    unsigned int test_partition = 0;
    for (unsigned int m = 0; m != 10; ++m) 
    {
        if (nodes[m].handle() != m)     ITFAILS;

        const Point_t &test_pt = nodes[m].coupling_object();
        if (test_pt.x() != counter)     ITFAILS;
        if (test_pt.y() != counter)     ITFAILS;
        if (test_pt.z() != counter)     ITFAILS;

        if (nodes[m].partition() != test_partition)   ITFAILS;
        if (nodes[m].get("power") != counter)         ITFAILS;
        if (nodes[m].get("temperature") != 2*counter) ITFAILS;
        if (nodes[m].get("volume") != 3*counter)      ITFAILS;

        counter += 1.0;
        test_partition += 2;
    }

    // Finally, print the info in a map node just to test the print function
    {
        // Make a point
        Point_t point(1.0, 2.5, 3.2);

        // Make a map node
        Map_Node_t node(4, 3, point);

        // Add some data
        node.add("power", 4.3);
        node.add("temperature", 2.1);

        // print the node
        node.print(cout);
        std::cout << std::endl;
    }

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "One processor Map_Node test ok on " << node;
        ut.passes( m.str() );
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Parallel_Unit_Test ut(argc, argv, denovo::release);

    node  = nemesis::node();
    nodes = nemesis::nodes();
    
    try
    {
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;

        if (nodes == 1)
        {
            serial_node_point_test(ut);

            gpass += ut.numPasses;
            gfail += ut.numFails;
            ut.reset();
        }

        else 
        {
            ut.passes("Scalar test only");

            gpass += ut.numPasses;
            gfail += ut.numFails;
            ut.reset();
        }
        
        // add up global passes and fails
        nemesis::global_sum(gpass);
        nemesis::global_sum(gfail);
        ut.numPasses = gpass;
        ut.numFails  = gfail;
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstMap_Node, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstMap_Node, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstMap_Node.cc
//---------------------------------------------------------------------------//
