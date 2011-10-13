//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/test/tstSTAR_Data.cc
 * \author Stuart R. Slattery
 * \date   Tue Jun 07 11:32:02 2011
 * \brief  Test for the Star_Data container
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
#include "utils/Packing_Utils.hh"
#include "release/Release.hh"
#include "mesh_type/Point.hh"
#include "../STAR_Data.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// make a vector of Point objects
std::vector< denovo::Point<int, double> > make_points()
{
    std::vector< denovo::Point<int, double> > points;

    double counter = 0.0;
    for (int i = 0; i < 10; ++i)
    {
        denovo::Point<int, double> new_point(counter, counter, counter, i);
        points.push_back(new_point);
        counter += 1.0;
    }

    return points;
}

//---------------------------------------------------------------------------//
// make a vector of volumes 
std::vector<double> make_volumes()
{
    std::vector<double> volumes;

    double counter = 0.0;
    for (int i = 0; i < 10; ++i)
    {
        volumes.push_back(counter);
        counter += 1.0;
    }

    return volumes;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void star_data_test(Parallel_Unit_Test &ut)
{
    // build the star data on node 0 and send it 
    if (node == 0)
    {
        // make a vector of points
        std::vector< denovo::Point<int, double> > points = make_points();
    
        // make a vector of volumes
        std::vector<double> volumes = make_volumes();
        UNIT_TEST(points.size() == volumes.size());

        // make a vector of Star_Data containers
        std::vector< coupler::STAR_Data<double> > star_data_vec;
        for (int i = 0; i < points.size(); ++i)
        {
            coupler::STAR_Data<double> new_star_data(points[i], volumes[i]);
            star_data_vec.push_back(new_star_data);
        }

        // create a packer
        denovo::Packer p;

        // compute the size of the message buffer
        p.compute_buffer_size_mode();
        for (int k = 0; k < star_data_vec.size(); ++k)
        {
            p << star_data_vec[k];
        }
        int buffer_size = p.size();
        UNIT_TEST(buffer_size > 0);

        // send a message to nodes 1, 2, and 3 to tell them how big the
        // incoming buffer is 
        for (int n = 0; n < nemesis::nodes(); ++n)
        {
            nemesis::send_async(&buffer_size, 1, n);
        }

        // create the buffer
        std::vector<char> buffer(buffer_size);
        p.set_buffer(buffer_size, &buffer[0]);

        // pack the star data
        for (int l = 0; l < star_data_vec.size(); ++l)
        {
            p << star_data_vec[l];
        }

        // send the buffer to nodes
        for (int n = 0; n < nemesis::nodes(); ++n)
        {
            nemesis::send_async(&buffer[0], buffer_size, n);
        }
    }

    // nodes get star data from node 0 and checks it
    // get a message from node 0 to see how big to make the buffer
    int incoming_size;
    nemesis::receive(&incoming_size, 1, 0);
        
    // setup a buffer to fill
    std::vector<char> incoming_buffer(incoming_size);
        
    // get the buffer from node 0
    nemesis::receive( &incoming_buffer[0], incoming_size, 0);

    // unpack the buffer into star data containers and check them
    double counter = 0.0;
    denovo::Unpacker u;
    u.set_buffer(incoming_size, &incoming_buffer[0]);
    for (int j = 0; j < 10; ++j)
    {
        denovo::Point<int, double> dummy_point(0.0, 0.0, 0.0, 0);
        coupler::STAR_Data<double> star_data(dummy_point, 0.0);

        u >> star_data;

        if ( !soft_equiv(star_data.point().x(), counter) )  ITFAILS;
        if ( !soft_equiv(star_data.point().y(), counter) )  ITFAILS;
        if ( !soft_equiv(star_data.point().z(), counter) )  ITFAILS;
        if (star_data.point().handle() != j)                ITFAILS;
        if ( !soft_equiv(star_data.volume(), counter) )     ITFAILS;

        counter += 1.0;
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Star_Data test ok on " << node;
        ut.passes(m.str());
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
        
        
        if(nodes > 1)
        {
            star_data_test(ut);
            gpass += ut.numPasses;
            gfail += ut.numFails;
            ut.reset();
        }
        else
        {
            ++gpass;
        }

        // add up global passes and fails
        nemesis::global_sum(gpass);
        nemesis::global_sum(gfail);
        ut.numPasses = gpass;
        ut.numFails  = gfail;
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing tstStar_Data, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstStar_Data, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstStar_Data.cc
//---------------------------------------------------------------------------//
