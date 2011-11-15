//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/test/tstMessenger.cc
 * \author Stuart Slattery
 * \date   Thu Jun 02 09:10:58 2011
 * \brief  Tests the Messenger class.
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
#include "utils/Packing_Utils.hh"
#include "../Messenger.hh"

using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

using denovo::SP;
using denovo::Packer;
using denovo::Unpacker;
using denovo::Point;

using database::Std_DB;
using kba::Simple_Partitioner;
using kba::Partitioner;

using coupler::Map_Node;
using coupler::Point_Map;
using coupler::Messenger;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
// Figure out the spatial partitioning
std::pair<int, int> make_spatial_partitioning(Parallel_Unit_Test &ut)
{
    UNIT_TEST(nemesis::nodes() == 1 || nemesis::nodes() == 2 ||
              nemesis::nodes() == 4);

    int num_blocks_i, num_blocks_j;

    if(nemesis::nodes() == 1)
    {
        num_blocks_i = 1;
        num_blocks_j = 1;
    }
    else if(nemesis::nodes() == 2)
    {
        num_blocks_i = 2;
        num_blocks_j = 1;
    }
    else if(nemesis::nodes() == 4)
    {
        num_blocks_i = 2;
        num_blocks_j = 2;
    }
    else
    {
        std::cout << "Wrong number of nodes: " << nemesis::nodes() << std::endl;
        UNIT_TEST(false);
    }

    return std::make_pair(num_blocks_i, num_blocks_j);
}

//---------------------------------------------------------------------------//
// Make a partitioner
SP<Partitioner> make_partitioner(int num_blocks_i, int num_blocks_j,
                                 Parallel_Unit_Test &ut)
{
    typedef kba::Simple_Partitioner             Simple_Partitioner;
    typedef Simple_Partitioner::SP_Std_DB       SP_Std_DB;
    typedef database::Std_DB                    Std_DB;
    typedef SP<Simple_Partitioner>              SP_Partitioner;

    // database
    Simple_Partitioner::SP_Std_DB db( new Std_DB("partition") );
    UNIT_TEST(db);

    // set data
    {
        db->new_key("num_blocks_i", num_blocks_i);
        db->new_key("num_blocks_j", num_blocks_j);
        db->new_key("num_cells_i", 10);
        db->new_key("num_cells_j", 10);
        db->new_key("num_cells_k", 10);
        db->new_key("delta_x", 1.0);
        db->new_key("delta_y", 1.0);
        db->new_key("delta_z", 1.0);
        db->new_key("num_z_blocks", 1);
        db->new_key("num_groups", 1);
    }

    // Instantiate the partitioner
    SP_Partitioner partitioner( new Simple_Partitioner(db) );
    UNIT_TEST(partitioner);

    // Do the build
    partitioner->build();

    return partitioner;
}

//---------------------------------------------------------------------------//
// Build a Point_Map
template<class OrdinateType_T, class DataType_T>
SP< Point_Map<OrdinateType_T, DataType_T> >
make_point_map(SP<Partitioner> partitioner, Parallel_Unit_Test &ut)
{
    typedef Point<OrdinateType_T, DataType_T>       Point_t;
    typedef Point_Map<OrdinateType_T, DataType_T>   Point_Map_t;
    typedef SP<Point_Map_t>                         SP_Point_Map;

    UNIT_TEST (partitioner);

    // Create a point map
    SP<Point_Map_t> point_map( new Point_Map_t(partitioner) );
    UNIT_TEST (point_map);

    // Make a bunch of points and add them to the map
    DataType_T x_inc = 0.5;
    DataType_T y_inc = 0.5;
    DataType_T z_inc = 0.5;
    DataType_T dx = 1.0;
    DataType_T dy = 1.0;
    DataType_T dz = 1.0;

    OrdinateType_T handle = 0;
    for(OrdinateType_T x = 0; x < 10; ++x)
    {
        for(OrdinateType_T y = 0; y < 10; ++y)
        {
            for(OrdinateType_T z = 0; z < 10; ++z)
            {
                DataType_T x_loc = dx * x + x_inc;
                DataType_T y_loc = dy * y + y_inc;
                DataType_T z_loc = dz * z + z_inc;
            
                // Make the point
                Point_t point(x_loc, y_loc, z_loc, handle);
                ++handle;
                
                // Add the point to the map
                point_map->add_point(point);
            }
        }
    }

    // Complete the mapping
    point_map->complete();

    return point_map;
}                


//---------------------------------------------------------------------------//
// Fill point map
template<class OrdinateType_T, class DataType_T>
void fill_point_map(SP< Point_Map<OrdinateType_T, DataType_T> > point_map,
                    Parallel_Unit_Test &ut)
{
    typedef Point_Map<OrdinateType_T, DataType_T>       Point_Map_t;
    typedef typename Point_Map_t::Map_Node_t            Map_Node_t;

    UNIT_TEST(point_map);

    for(OrdinateType_T index = 0; index < point_map->num_nodes(); ++index)
    {
        // Get the node
        Map_Node_t& node = point_map->get_node(index);

        // If we own this node, insert a power (equal to the index)
        node.register_data("power", OrdinateType_T(index) );

        // Verify that the power was inserted correctly
        UNIT_TEST( soft_equiv( node.get_data("power"), DataType_T(index) ) );
    }
}

//---------------------------------------------------------------------------//
// Check point map
template<class OrdinateType_T, class DataType_T>
void check_point_map(SP< Point_Map<OrdinateType_T, DataType_T> > point_map,
                     Parallel_Unit_Test &ut)
{
    typedef Point_Map<OrdinateType_T, DataType_T>       Point_Map_t;
    typedef typename Point_Map_t::Map_Node_t            Map_Node_t;

    UNIT_TEST(point_map);

    for(OrdinateType_T index = 0; index < point_map->num_nodes(); ++index)
    {
        // Get the node
        Map_Node_t& node = point_map->get_node(index);

        // Check that we have a power (every node should be filled now)
        UNIT_TEST( node.data_exists("power") );

        // Check that this node contains a power equal to the index
        UNIT_TEST( node.get_data("power") == DataType_T(index) );
    }
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

// messenger test
void messenger_test(Parallel_Unit_Test &ut)
{
    typedef int                                 OrdinateType;
    typedef double                              DataType;
    typedef SP<kba::Partitioner>                SP_Partitioner;
    typedef Point_Map<OrdinateType, DataType>   Point_Map_t;
    typedef SP<Point_Map_t>                     SP_Point_Map;

    // Make the spatial partitioning
    std::pair<int, int> spatial_partitioning = make_spatial_partitioning(ut);

    // Make the partitioner
    SP_Partitioner partitioner = 
        make_partitioner(spatial_partitioning.first, 
                         spatial_partitioning.second, ut);
    UNIT_TEST(partitioner);

    // Make a point map
    SP_Point_Map point_map = 
        make_point_map<OrdinateType, DataType>(partitioner, ut);
    UNIT_TEST(point_map);

    // Fill point map
    fill_point_map<OrdinateType, DataType>(point_map, ut);

    // Get the global communicator
    nemesis::Communicator_t global_comm = get_comm_world();

    // Make the messenger
    Messenger<DataType> messenger(global_comm, point_map);

    // Do the communication
    messenger.communicate("power");

    // Check the data
    check_point_map<OrdinateType, DataType>(point_map, ut);

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "Messenger test ok on " << nemesis::node();
        ut.passes( m.str() );
    }

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

        if(nemesis::nodes() > 1)
        {
            messenger_test(ut);
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
        std::cout << "ERROR: While testing tstMessenger, " 
                  << err.what()
                  << std::endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstMessenger, " 
                  << "An unknown exception was thrown."
                  << std::endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstMessenger.cc
//---------------------------------------------------------------------------//
