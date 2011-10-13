//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/test/tstStar_Denovo_Coupling.cc
 * \author Stuart R. Slattery
 * \date   Thu Jun 09 16:47:01 2011
 * \brief  Test for Star-Denovo coupling
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template_c4_test.cc,v 1.7 2008/01/02 22:50:26 9te Exp $
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
#include <limits>

#include "harness/DBC.hh"
#include "harness/Soft_Equivalence.hh"
#include "comm/global.hh"
#include "comm/Parallel_Unit_Test.hh"
#include "release/Release.hh"
#include "utils/SP.hh"
#include "neutronics/Neutronics.hh"
#include "../Coupler.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

using coupler::Coupler;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
nemesis::Communicator_t get_comm_world()
{
#ifdef COMM_MPI
    return MPI_COMM_WORLD;
#else
    return 1;
#endif
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// one processor unit test for Star-Denovo coupling
void coupler_test(int num_I_blocks, int num_J_blocks, int num_sets, 
                  Parallel_Unit_Test &ut)
{
    // typedefs
    typedef denovo::SP<neutronics::Neutronics>          SP_Neutronics;
    typedef neutronics::Neutronics::const_View_Field    const_View_Field;

    // Create the denovo input filename
    std::ostringstream m;
    m << "NeutronicsTest_" << num_I_blocks << "_" << num_J_blocks << "_"
      << num_sets << ".in";
    std::string filename = m.str();

    // make a Neutronics object
    SP_Neutronics neutronics( new neutronics::Neutronics("dummy", filename) );
    UNIT_TEST(neutronics);

    // neutronics setup
    neutronics->setup();

    // make a Coupler on all nodes
    Coupler coupler(get_comm_world(), get_comm_world(), get_comm_world(), 
                    neutronics, neutronics);    

    // assign the neutronics object to the coupler
    coupler.register_neutronics(neutronics);

    // build the transfer map
    coupler.build_map("STAR_test_geom.inp");        

    // do transport
    neutronics->transport();

    // Get the intensities
    const_View_Field intensities = neutronics->get_intensities();

    // Print the intensities
    /*
    for(int i=0; i < intensities.size(); ++i)
    {
        std::cout << intensities[i] << std::endl;
    }
    */

    // write out the Star data file with powers
    coupler.transfer_power("STAR_power_file.out");

    // check the results
    // reference power value
    int num_groups = 4;
    double ref_power = num_groups * (205.0) * sqrt(4.0 * 3.141592654) * 
                       (100.0 * 100.0 * 100.0) / 
                       6.24150974e12;

    std::cout << "Ref power = " << ref_power << std::endl;

    // read the Star output file back in to check it
    if (node == 0)
    {
        vector<double> in_x_vec;
        vector<double> in_y_vec;
        vector<double> in_z_vec;
        
        vector<double> out_x_vec;
        vector<double> out_y_vec;
        vector<double> out_z_vec;
        vector<double> power_vec;

        std::ifstream input("STAR_test_geom.inp");
        if (!input) ITFAILS;

        std::ifstream output("STAR_power_file.out");
        if (!output) ITFAILS;

        // Skip the headers
        input.ignore(std::numeric_limits<int>::max(), '\n');
        output.ignore(std::numeric_limits<int>::max(), '\n');

        double dummy, volume, in_x, in_y, in_z;
        double power, out_x, out_y, out_z;

        while(!input.eof())
        {
            input >> dummy >> dummy >> dummy >> volume >> in_x >> in_y >> in_z;
            output >> power >> out_x >> out_y >> out_z;

            in_x_vec.push_back(in_x);
            in_y_vec.push_back(in_y);
            in_z_vec.push_back(in_z);
            out_x_vec.push_back(out_x);
            out_y_vec.push_back(out_y);
            out_z_vec.push_back(out_z);
            power_vec.push_back(power);
        }

        // Now, loop over the powers and check the vals
        for(int i = 0; i < power_vec.size(); ++i)
        {
            // Check the power against the reference
            if ( !soft_equiv(power, ref_power, 1.0e-5) ) ITFAILS;

            // Loop over the input and see if we can find the output x
            int j = 0;
            for(; j < in_x_vec.size(); ++j)
            {
                if( soft_equiv(in_x_vec[j], out_x_vec[i]) )
                    break;
            }

            if(j == in_x_vec.size())     ITFAILS;

            if( !soft_equiv(in_x_vec[j], out_x_vec[i]) )    ITFAILS;
            if( !soft_equiv(in_y_vec[j], out_y_vec[i]) )    ITFAILS;
            if( !soft_equiv(in_z_vec[j], out_z_vec[i]) )    ITFAILS;
        }
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Star_Denovo_Coupling test with " << num_I_blocks << " I blocks, "
          << num_J_blocks << " J blocks, and " << num_sets << " sets ok on " 
          << node;
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

        if (nodes == 1)
        {
            /*
            int num_I_blocks = 1;
            int num_J_blocks = 1;
            int num_sets = 1;
            coupler_test(num_I_blocks, num_J_blocks, num_sets, ut);
            */
            ++gpass;
        }
        else if(nodes == 2)
        {
            int num_I_blocks = 2;
            int num_J_blocks = 1;
            int num_sets = 1;
            coupler_test(num_I_blocks, num_J_blocks, num_sets, ut);
            
            num_I_blocks = 1;
            num_J_blocks = 1;
            num_sets = 2;
            coupler_test(num_I_blocks, num_J_blocks, num_sets, ut);
        }
        else if(nodes == 4)
        {
            int num_I_blocks = 2;
            int num_J_blocks = 2;
            int num_sets = 1;
            coupler_test(num_I_blocks, num_J_blocks, num_sets, ut);
            
            num_I_blocks = 2;
            num_J_blocks = 1;
            num_sets = 2;
            coupler_test(num_I_blocks, num_J_blocks, num_sets, ut);

            num_I_blocks = 4;
            num_J_blocks = 1;
            num_sets = 1;
            coupler_test(num_I_blocks, num_J_blocks, num_sets, ut);
            
            num_I_blocks = 1;
            num_J_blocks = 1;
            num_sets = 4;
            coupler_test(num_I_blocks, num_J_blocks, num_sets, ut);
        }
        else
        {
            ++gpass;
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
        std::cout << "ERROR: While testing tstStar_Denovo_Coupling, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstStar_Denovo_Coupling, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstStar_Denovo_Coupling.cc
//---------------------------------------------------------------------------//
