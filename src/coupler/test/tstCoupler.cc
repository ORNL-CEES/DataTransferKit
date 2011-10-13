//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/test/tstCoupler.cc
 * \author Stuart Slattery
 * \date   Wed Jun 15 09:15:03 2011
 * \brief  Unit tests for Coupler
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
#include "utils/SP.hh"
#include "mesh_type/Point.hh"
#include "neutronics/Neutronics.hh"
#include "../Coupler.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

using neutronics::Neutronics;

using coupler::Coupler;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// pseudo CFD application to hold quadrature points and thermal sources
class CFD_app
{
  public:

    CFD_app() { /* ... */ }

    ~CFD_app() { /* ... */ }

    void set_handles(vector<int> handles) { d_handles = handles; }

    vector<int> get_handles() { return d_handles; }

    void set_points(vector<double> points) { d_points = points; }

    vector<double> get_points() { return d_points; }

    void set_qqq(vector<double> qqq) { d_qqq = qqq; }

    vector<double> get_qqq() { return d_qqq; }

  private:

    vector<int> d_handles;

    vector<double> d_points;

    vector<double> d_qqq;
};

//---------------------------------------------------------------------------//
// make a set of interleaved point coordinates
vector<double> make_points()
{
    vector<double> points(12, 0.0);

    // point 1
    points[0] = 0.0043;
    points[1] = 1.322;
    points[2] = 5.3365;

    // point 2
    points[3] = 9.56;
    points[4] = 7.56;
    points[5] = 2.36;

    // point 3
    points[6] = 2.654;
    points[7] = 8.225;
    points[8] = 5.5556;
    
    // point 4
    points[9] = 6.5458;
    points[10] = 3.748;
    points[11] = 7.1645;

    return points;
}

//---------------------------------------------------------------------------//
// make a set of handles for the points we created
vector<int> make_handles()
{
    vector<int> handles(4, 0);

    for (vector<int>::iterator it = handles.begin(),
                           it_end = handles.end();
         it != it_end; ++it)
    {
        *it = it - handles.begin();
    }

    return handles;
}

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
// test coupler data transfer with both applications on every node
void transfer_test(int num_I_blocks, int num_J_blocks, int num_sets,
                   Parallel_Unit_Test &ut)
{
    // typedefs
    typedef denovo::SP<CFD_app>                         SP_CFD;
    typedef denovo::SP<Neutronics>                      SP_Neutronics;
    typedef nemesis::Communicator_t                     Communicator_t;


    // NEUTRONICS SETUP

    // Set the neutronics communicator
    nemesis::set_internal_comm( get_comm_world() );
    Communicator_t comm_neutronics;
    nemesis::split(0, nemesis::nodes() - nemesis::node() - 1, comm_neutronics);
    nemesis::set_default(comm_neutronics);

    // Create the denovo input filename
    std::ostringstream m;
    m << "NeutronicsTest_" << num_I_blocks << "_" << num_J_blocks << "_"
      << num_sets << ".in";
    std::string filename = m.str();

    // make a Neutronics object
    SP_Neutronics neutronics( new Neutronics("dummy", filename) );
    UNIT_TEST(neutronics);

    // neutronics setup
    neutronics->setup();

    // barrier after neutronics setup
    nemesis::global_barrier();

    // CFD SETUP

    // create a cfd communicator
    nemesis::set_internal_comm( get_comm_world() );
    Communicator_t comm_cfd;
    nemesis::split(0, nemesis::node(), comm_cfd);
    nemesis::set_internal_comm(comm_cfd);

    // make a CFD object
    SP_CFD cfd( new CFD_app() );
    UNIT_TEST(cfd);

    // assign some quadrature points and handles to the cfd application
    vector<double> points = make_points();
    cfd->set_points(points);
    vector<int> handles = make_handles();
    cfd->set_handles(handles);

    // barrier after cfd setup
    nemesis::global_barrier();
    nemesis::set_internal_comm( get_comm_world() );


    // COUPLER SETUP

    // setup the coupler
    Coupler coupler(get_comm_world(), comm_neutronics, comm_cfd, 
                    neutronics, cfd);

    // neutronics application nodes register a neutronics object 
    if (neutronics) 
    {
        coupler.register_neutronics(neutronics);
    }

    // assign points from CFD to the coupler
    if (cfd)
    {
        coupler.register_points( cfd->get_points(), cfd->get_handles() );
    }

    // build the transfer map
    coupler.build_map();

    // do transport
    if (neutronics) 
    {
        neutronics->transport();
    }

    // transfer powers from neutronics to cfd
    coupler.transfer_power();

    // assign the transfered powers to cfd and check the results
    if (cfd) 
    {
        vector<double> powers(cfd->get_handles().size(), 0.0);

        coupler.get_power( powers );

        cfd->set_qqq( powers );

        vector<double> cfd_powers = cfd->get_qqq();

        int num_groups = 4;
        double ref_power = num_groups * (205.0) * sqrt(4.0 * 3.141592654) * 
                           (100.0 * 100.0 * 100.0) / 
                           6.24150974e12;

        for (vector<double>::const_iterator it = cfd_powers.begin(),
                                        it_end = cfd_powers.end();
             it != it_end; ++it)
        {
            cout << "Power: " << *it << " Ref: " << ref_power << endl;
            UNIT_TEST( soft_equiv( *it, ref_power, 1.0e-5 ) );
        }
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Coupler transfer test with " << num_I_blocks << " I blocks, "
          << num_J_blocks << " J blocks, and " << num_sets << " sets ok on " 
          << node;
        ut.passes(m.str());
    }

    nemesis::global_barrier();    
    nemesis::set_default( get_comm_world() );
}

//---------------------------------------------------------------------------//
// Unit test for varying communicator sizes. Here the applications are split
// across the parallel domain. E.g. for a 4 process problem, 1 processor holds
// neutronics while the other three hold a cfd application.
void variable_comm_test(int num_I_blocks, int num_J_blocks, int num_sets,
                        int num_neutronics, int num_cfd, Parallel_Unit_Test &ut)
{
    // typedefs
    typedef denovo::SP<CFD_app>                         SP_CFD;
    typedef denovo::SP<Neutronics>                      SP_Neutronics;
    typedef nemesis::Communicator_t                     Communicator_t;

    // initialize application pointers
    SP_Neutronics neutronics;
    SP_CFD cfd;

    int global_id = nemesis::node();

    // setup communicators
    Communicator_t comm_split;
    if (global_id < num_neutronics)
    {
        nemesis::split(0, 0, comm_split);
    }
    else
    {
        nemesis::split(1, 0, comm_split);
    }

    std::cout << "Ok_a" << std::endl;

    Communicator_t comm_neutronics;
#ifdef COMM_MPI
    int result = MPI_Comm_dup(comm_split, &comm_neutronics);
    UNIT_TEST(result == MPI_SUCCESS);
    std::cout << "Ok_a.5" << std::endl;
#endif

    Communicator_t comm_cfd;
#ifdef COMM_MPI
    result = MPI_Comm_dup(comm_split, &comm_cfd);
    UNIT_TEST(result == MPI_SUCCESS);
#endif

    std::cout << "Ok_b" << std::endl;

    // NEUTRONICS SETUP
    nemesis::set_default(comm_neutronics);

    std::cout << "Ok_c" << std::endl;

    if (global_id < num_neutronics)
    {
        std::ostringstream m;
        m << "NeutronicsTest_" << num_I_blocks << "_" << num_J_blocks << "_"
          << num_sets << ".in";
        std::string filename = m.str();

        neutronics = new Neutronics("dummy", filename);
        UNIT_TEST(neutronics);
        neutronics->setup();
    }

    // CFD SETUP
    else
    {
        nemesis::set_internal_comm( get_comm_world() );

        // make a CFD object
        cfd = new CFD_app();
        UNIT_TEST(cfd);

        // assign some quadrature points and handles to the cfd application
        vector<double> points = make_points();
        cfd->set_points(points);
        vector<int> handles = make_handles();
        cfd->set_handles(handles);

        nemesis::reset_internal_comm();
    }

    std::cout << "Ok_d" << std::endl;


    // barrier after setup
    nemesis::global_barrier();
    nemesis::set_internal_comm( get_comm_world() );

    // COUPLER SETUP

    // setup the coupler
    Coupler coupler(get_comm_world(), comm_split, comm_split, 
                    neutronics, cfd);

    // neutronics application nodes register a neutronics object 
    if (neutronics) 
    {
        coupler.register_neutronics(neutronics);
    }

    // assign points from CFD to the coupler
    if (cfd)
    {
        coupler.register_points( cfd->get_points(), cfd->get_handles() );
    }

    // build the transfer map
    coupler.build_map();

    // do transport
    if (neutronics) 
    {
        neutronics->transport();
    }

    nemesis::global_barrier();
    std::cout << "Ok_e" << std::endl;

    // transfer powers from neutronics to cfd
    coupler.transfer_power();

    // assign the transfered powers to cfd and check the results
    if (cfd) 
    {
        vector<double> cfd_powers(cfd->get_handles().size(), 0.0);
        coupler.get_power( cfd_powers );

        int num_groups = 4;
        double ref_power = num_groups * (205.0) * sqrt(4.0 * 3.141592654) * 
                           (100.0 * 100.0 * 100.0) / 
                           6.24150974e12;

        for (vector<double>::const_iterator it = cfd_powers.begin(),
                                        it_end = cfd_powers.end();
             it != it_end; ++it)
        {
            UNIT_TEST( soft_equiv( *it, ref_power, 1.0e-5 ) );
            cout << "Transferred: " << *it << " Ref: " << ref_power << endl;
        }
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Coupler variable comm test with " << num_I_blocks << " I blocks, "
          << num_J_blocks << " J blocks, and " << num_sets << " sets ok on " 
          << node;
        ut.passes(m.str());
    }

    nemesis::global_barrier();    
    nemesis::set_default( get_comm_world() );
}

//---------------------------------------------------------------------------//
// Unit test for overlapping communicators.
void overlap_comm_test(int num_I_blocks, int num_J_blocks, int num_sets,
                       int num_neutronics, int num_cfd, Parallel_Unit_Test &ut)
{
    // this test is for four processors only
    if (nemesis::nodes() != 4) 
        return;

    // typedefs
    typedef denovo::SP<CFD_app>                         SP_CFD;
    typedef denovo::SP<Neutronics>                      SP_Neutronics;
    typedef nemesis::Communicator_t                     Communicator_t;

    // initialize application pointers
    SP_Neutronics neutronics;
    SP_CFD cfd;

    int global_id = nemesis::node();

    // setup communicators
    Communicator_t comm_split;
    if (nemesis::node() < 2)
    {
        nemesis::split(0, 0, comm_split);
    }
    else
    {
        nemesis::split(1, 0, comm_split);
    }

    Communicator_t comm_neutronics = get_comm_world();
    /*
      int result = MPI_Comm_dup(get_comm_world(), &comm_neutronics);
      UNIT_TEST(result == MPI_SUCCESS);
    */

    Communicator_t comm_cfd = comm_split;
    /*
      result = MPI_Comm_dup(comm_split, &comm_cfd);
      UNIT_TEST(result == MPI_SUCCESS);
    */

    // NEUTRONICS SETUP
    nemesis::set_default(comm_neutronics);
    nemesis::set_internal_comm( get_comm_world() );

    std::ostringstream m;
    m << "NeutronicsTest_" << num_I_blocks << "_" << num_J_blocks << "_"
      << num_sets << ".in";
    std::string filename = m.str();

    neutronics = new Neutronics("dummy", filename);
    UNIT_TEST(neutronics);
    neutronics->setup();

    // CFD SETUP
    if (nemesis::node() < 2)
    {
        // make a CFD object
        cfd = new CFD_app();
        UNIT_TEST(cfd);

        // assign some quadrature points and handles to the cfd application
        vector<double> points = make_points();
        cfd->set_points(points);
        vector<int> handles = make_handles();
        cfd->set_handles(handles);
    }

    // barrier after setup
    nemesis::global_barrier();

    // COUPLER SETUP

    // setup the coupler
    Coupler coupler(get_comm_world(), comm_neutronics, comm_cfd, 
                    neutronics, cfd);
 
    // neutronics application nodes register a neutronics object 
    if (neutronics) 
    {
        coupler.register_neutronics(neutronics);
    }

    // assign points from CFD to the coupler
    if (cfd)
    {
        coupler.register_points( cfd->get_points(), cfd->get_handles() );
    }

    // build the transfer map
    coupler.build_map();

    // do transport
    if (neutronics) 
    {
        neutronics->transport();
    }

    // transfer powers from neutronics to cfd
    coupler.transfer_power();

    // assign the transfered powers to cfd and check the results
    if (cfd) 
    {
        vector<double> cfd_powers(cfd->get_handles().size(), 0.0);
        coupler.get_power( cfd_powers );

        int num_groups = 4;
        double ref_power = num_groups * (205.0) * sqrt(4.0 * 3.141592654) * 
                           (100.0 * 100.0 * 100.0) / 
                           6.24150974e12;

        for (vector<double>::const_iterator it = cfd_powers.begin(),
                                        it_end = cfd_powers.end();
             it != it_end; ++it)
        {
            UNIT_TEST( soft_equiv( *it, ref_power, 1.0e-5 ) );
            cout << "Transferred: " << *it << " Ref: " << ref_power << endl;
        }
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Coupler variable comm test with " << num_I_blocks << " I blocks, "
          << num_J_blocks << " J blocks, and " << num_sets << " sets ok on " 
          << node;
        ut.passes(m.str());
    }

    nemesis::global_barrier();
    nemesis::set_default( get_comm_world() );
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[])
{
    Parallel_Unit_Test ut(argc, argv, denovo::release);

    try
    {
        // >>> UNIT TESTS
        int gpass = 0;
        int gfail = 0;

        if (nemesis::nodes() == 1)
        {
            /*
            int num_I_blocks = 1;
            int num_J_blocks = 1;
            int num_sets = 1;
            transfer_test(num_I_blocks, num_J_blocks, num_sets, ut);
            */

            ++gpass;
        }
        else if(nemesis::nodes() == 2)
        {
            int num_I_blocks = 2;
            int num_J_blocks = 1;
            int num_sets = 1;
            //transfer_test(num_I_blocks, num_J_blocks, num_sets, ut);
            /*
            num_I_blocks = 1;
            num_J_blocks = 1;
            num_sets = 2;
            transfer_test(num_I_blocks, num_J_blocks, num_sets, ut);
            */
            num_I_blocks = 1;
            num_J_blocks = 1;
            num_sets = 1;
            int num_neutronics = 1;
            int num_cfd = 1;
            /*
            variable_comm_test(num_I_blocks, num_J_blocks, num_sets,
                               num_neutronics, num_cfd, ut);
            */
            
            ++gpass;
        }
        else if(nemesis::nodes() == 4)
        {
            int num_I_blocks = 2;
            int num_J_blocks = 2;
            int num_sets = 1;
            //transfer_test(num_I_blocks, num_J_blocks, num_sets, ut);
            /*
            num_I_blocks = 2;
            num_J_blocks = 1;
            num_sets = 2;
            transfer_test(num_I_blocks, num_J_blocks, num_sets, ut);
            */
            num_I_blocks = 4;
            num_J_blocks = 1;
            num_sets = 1;
            //transfer_test(num_I_blocks, num_J_blocks, num_sets, ut);
            /*            
            num_I_blocks = 1;
            num_J_blocks = 1;
            num_sets = 4;
            transfer_test(num_I_blocks, num_J_blocks, num_sets, ut);
            */
            num_I_blocks = 1;
            num_J_blocks = 1;
            num_sets = 1;
            int num_neutronics = 1;
            int num_cfd = 3;
            /*
            variable_comm_test(num_I_blocks, num_J_blocks, num_sets,
                               num_neutronics, num_cfd, ut);
            */

            num_I_blocks = 2;
            num_J_blocks = 1;
            num_sets = 1;
            num_neutronics = 2;
            num_cfd = 2;
            /*
            variable_comm_test(num_I_blocks, num_J_blocks, num_sets,
                               num_neutronics, num_cfd, ut);
            */

            num_I_blocks = 3;
            num_J_blocks = 1;
            num_sets = 1;
            num_neutronics = 3;
            num_cfd = 1;
            /*
            variable_comm_test(num_I_blocks, num_J_blocks, num_sets,
                               num_neutronics, num_cfd, ut);
            */

            num_I_blocks = 2;
            num_J_blocks = 2;
            num_sets = 1;
            num_neutronics = 4;
            num_cfd = 2;
            /*
            overlap_comm_test(num_I_blocks, num_J_blocks, num_sets,
                              num_neutronics, num_cfd, ut);
            */

            num_I_blocks = 4;
            num_J_blocks = 1;
            num_sets = 1;
            num_neutronics = 4;
            num_cfd = 2;
            /*
            overlap_comm_test(num_I_blocks, num_J_blocks, num_sets,
                              num_neutronics, num_cfd, ut);
            */
            ++gpass;
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
        std::cout << "ERROR: While testing tstCoupler, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstCoupler, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstCoupler.cc
//---------------------------------------------------------------------------//
