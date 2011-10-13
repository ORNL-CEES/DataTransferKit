//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/test/tstPoint_Map.cc
 * \author Stuart R. Slattery
 * \date   Thu May 26 21:52:37 2011
 * \brief  Tests for the Point_Map class
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
#include "utils/SP.hh"
#include "database/Std_DB.hh"
#include "release/Release.hh"
#include "kba_mesh/Partitioner.hh"
#include "kba_mesh/Simple_Partitioner.hh"
#include "../Point_Map.hh"

using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

using denovo::SP;

using database::Std_DB;
using kba::Simple_Partitioner;
using kba::Partitioner;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);


//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// create an interleaved set of field points
template<class DataType>
std::vector<DataType> make_field_pts()
{
    std::vector<DataType> pts;

    // first point
    pts.push_back(2.5);
    pts.push_back(2.5);
    pts.push_back(2.5);

    // second point
    pts.push_back(7.5);
    pts.push_back(2.5);
    pts.push_back(2.5);

    // third point
    pts.push_back(2.5);
    pts.push_back(7.5);
    pts.push_back(2.5);

    // fourth point
    pts.push_back(7.5);
    pts.push_back(7.5);
    pts.push_back(2.5);

    // fifth point (outside of global mesh)
    pts.push_back(-2.5);
    pts.push_back(2.5);
    pts.push_back(2.5);

    return pts;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

// serial test
void one_pe_map_test(Parallel_Unit_Test &ut)
{
    typedef int                                         OrdinateType;
    typedef double                                      DataType;
    typedef coupler::Point_Map<OrdinateType, DataType>  Point_Map_t;
    typedef typename Point_Map_t::Point_t               Point_t;

    // this test is serial
    if (nodes > 1)
        return;

    // make a partitioner
    SP<Partitioner> p;

    // database
    Simple_Partitioner::SP_Std_DB db( new Std_DB("partition") );
    UNIT_TEST(db);

    // set data
    {
        db->new_key("num_blocks_i", 1);
        db->new_key("num_blocks_j", 1);
        db->new_key("num_cells_i", 10);
        db->new_key("num_cells_j", 10);
        db->new_key("num_cells_k", 10);
        db->new_key("delta_x", 1.0);
        db->new_key("delta_y", 1.0);
        db->new_key("delta_z", 1.0);
        db->new_key("num_z_blocks", 1);
        db->new_key("num_groups", 1);
    }

    // make the simple partitioner
    {
        // simple partitioner specialization
        p = new Simple_Partitioner(db);
        UNIT_TEST(p);

        if (!p) ITFAILS;
    }
    
    // build the mesh
    p->build();

    // pretend that we are on a CFD application node
    if (node == 0) 
    {
        // make a point map
        Point_Map_t map(p);

        // make a field of points
        std::vector<DataType> field_pts = make_field_pts<DataType>();
	
        // add the points field to the map
        for (OrdinateType i = 0; i < 5; ++i)
        { 
            map.add_point( Point_t(field_pts[3*i], field_pts[3*i+1], 
                                   field_pts[3*i+2], i) );
        }	

        // make sure all of the field points got into the map
        if (map.num_nodes() != 5) ITFAILS;
	
        // complete the map
        map.complete();
	
        // check the status of the map
        if (!map.status()) ITFAILS;
	
        // insert some fake power data into the map
        {
            // point 1
            map.insert_data(0, "power", 1.35);
	    
            // point 2
            map.insert_data(1, "power", 8.964);
	    
            // point 3
            map.insert_data(2, "power", 0.201);
	    
            // point 4
            map.insert_data(3, "power", 7.775);
	    
            // point 5
            map.insert_data(4, "power", 982.2);
        }
	
        // verify the contents of the map
        {
            // number of nodes in the map
            if (map.num_nodes() != 5)           ITFAILS;
	    
            // unique set of partition ids
            if (map.partitions().size() != 2)   ITFAILS;
            if (map.partitions()[0] != -1)      ITFAILS;
            if (map.partitions()[1] != 0)       ITFAILS;
	    
            // point 1
            if (map.get_node(0).partition() != 0)                   ITFAILS;
            if ( !soft_equiv(map.get_data(0, "power"), 1.35) )      ITFAILS;
            if ( !soft_equiv(map.get_node(0).point().x(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(0).point().y(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(0).point().z(), 2.5) )    ITFAILS;
            if( map.get_node(0).point().handle() != 0)              ITFAILS;
	    
            // point 2
            if (map.get_node(1).partition() != 0)                   ITFAILS;
            if ( !soft_equiv(map.get_data(1, "power"), 8.964) )     ITFAILS;
            if ( !soft_equiv(map.get_node(1).point().x(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(1).point().y(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(1).point().z(), 2.5) )    ITFAILS;
            if(map.get_node(1).point().handle() != 1)               ITFAILS;
	    
            // point 3
            if (map.get_node(2).partition() != 0 )                  ITFAILS;
            if ( !soft_equiv(map.get_data(2, "power"), 0.201) )     ITFAILS;
            if ( !soft_equiv(map.get_node(2).point().x(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(2).point().y(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(2).point().z(), 2.5) )    ITFAILS;
            if(map.get_node(2).point().handle() != 2)               ITFAILS;
	    
            // point 4
            if (map.get_node(3).partition() != 0 )                  ITFAILS;
            if ( !soft_equiv(map.get_data(3, "power"), 7.775) )     ITFAILS;

            if ( !soft_equiv(map.get_node(3).point().x(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(3).point().y(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(3).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(3).point().handle() != 3)              ITFAILS;

            // point 5
            if (map.get_node(4).partition() != -1 )                 ITFAILS;
            if ( !soft_equiv(map.get_data(4, "power"), 982.2) )     ITFAILS;

            if ( !soft_equiv(map.get_node(4).point().x(), -2.5) )   ITFAILS;
            if ( !soft_equiv(map.get_node(4).point().y(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(4).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(4).point().handle() != 4)              ITFAILS;
        }

        {
            // Try printing the map
            map.print(std::cout);
        }
    }

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "One processor Point_Map test ok on " << node;
        ut.passes( m.str() );
    }
}

//---------------------------------------------------------------------------//
// two processor test - 2 I blocks, 1 J block
void two_pe_map_test(Parallel_Unit_Test &ut)
{
    typedef int                                         OrdinateType;
    typedef double                                      DataType;
    typedef coupler::Point_Map<OrdinateType, DataType>  Point_Map_t;
    typedef typename Point_Map_t::Point_t               Point_t;

    // this test is for two nodes
    if (nodes != 2)
        return;

    // make a partitioner
    SP<Partitioner> p;

    // database
    Simple_Partitioner::SP_Std_DB db(new Std_DB("partition"));
    UNIT_TEST(db);

    // set data
    {
        db->new_key("num_blocks_i", 2);
        db->new_key("num_blocks_j", 1);
        db->new_key("num_cells_i", 10);
        db->new_key("num_cells_j", 10);
        db->new_key("num_cells_k", 10);
        db->new_key("delta_x", 1.0);
        db->new_key("delta_y", 1.0);
        db->new_key("delta_z", 1.0);
        db->new_key("num_z_blocks", 1);
        db->new_key("num_groups", 1);
    }

    // make the simple partitioner
    {
        // simple partitioner specialization
        p = new Simple_Partitioner(db);
        UNIT_TEST(p);

        if (!p) ITFAILS;
    }
    
    // build the mesh
    p->build();

    // pretend we are on a CFD application node
    if (node == 0) 
    {	
        // make a point map
        Point_Map_t map(p);
	
        // make a field of points
        std::vector<DataType> field_pts = make_field_pts<DataType>();
	
        // add the point field to the map
        for (OrdinateType i = 0; i < 5; ++i)
        {
            map.add_point( Point_t(field_pts[3*i], field_pts[3*i+1], 
                                   field_pts[3*i+2], i) );
        }

        // make sure all of the field points got into the map
        if (map.num_nodes() != 5) ITFAILS;
	
        // complete the map
        map.complete();
	
        // check the status of the map
        if (!map.status()) ITFAILS;
	
        // insert some fake data into the map
        {
            // point 1
            map.insert_data(0, "power", 1.35);
	    
            // point 2
            map.insert_data(1, "power", 8.964);
	    
            // point 3
            map.insert_data(2, "power", 0.201);
	    
            // point 4
            map.insert_data(3, "power", 7.775);
	    
            // point 5
            map.insert_data(4, "power", 982.2);
        }
	
        // verify the contents of the map
        {
            // number of nodes in the map
            if (map.num_nodes() != 5)               ITFAILS;
	    
            // unique set of partition ids
            if (map.partitions().size() != 3)       ITFAILS;
            if (map.partitions()[0] != -1)          ITFAILS;
            if (map.partitions()[1] != 0)           ITFAILS;
            if (map.partitions()[2] != 1)           ITFAILS;
	    
            // point 1
            if (map.get_node(0).partition() != 0)                   ITFAILS;
            if ( !soft_equiv(map.get_data(0, "power"), 1.35) )      ITFAILS;
            if ( !soft_equiv(map.get_node(0).point().x(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(0).point().y(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(0).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(0).point().handle() != 0)              ITFAILS;
	    
            // point 2
            if (map.get_node(1).partition() != 1)                   ITFAILS;
            if ( !soft_equiv(map.get_data(1, "power"), 8.964) )     ITFAILS;
            if ( !soft_equiv(map.get_node(1).point().x(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(1).point().y(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(1).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(1).point().handle() != 1)              ITFAILS;
    
            // point 3
            if (map.get_node(2).partition() != 0)                   ITFAILS;
            if ( !soft_equiv(map.get_data(2, "power"), 0.201) )     ITFAILS;
            if ( !soft_equiv(map.get_node(2).point().x(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(2).point().y(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(2).point().z(), 2.5) )    ITFAILS;
            if(map.get_node(2).point().handle() != 2)               ITFAILS;
	    
            // point 4
            if (map.get_node(3).partition() != 1)                   ITFAILS;
            if ( !soft_equiv(map.get_data(3, "power"), 7.775) )     ITFAILS;
            if ( !soft_equiv(map.get_node(3).point().x(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(3).point().y(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(3).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(3).point().handle() != 3)              ITFAILS;

            // point 5
            if (map.get_node(4).partition() != -1)                  ITFAILS;
            if ( !soft_equiv(map.get_data(4, "power"), 982.2) )     ITFAILS;
            if ( !soft_equiv(map.get_node(4).point().x(), -2.5) )   ITFAILS;
            if ( !soft_equiv(map.get_node(4).point().y(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(4).point().z(), 2.5) )    ITFAILS;
            if( map.get_node(4).point().handle() != 4)              ITFAILS;
        }
    }

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "Two processor Point_Map test ok on " << node;
        ut.passes( m.str() );
    }
}

//---------------------------------------------------------------------------//
// two processor test - 1 I block, 2 J blocks
void two_pe_switch_map_test(Parallel_Unit_Test &ut)
{
    typedef int                                         OrdinateType;
    typedef double                                      DataType;
    typedef coupler::Point_Map<OrdinateType, DataType>  Point_Map_t;
    typedef typename Point_Map_t::Point_t               Point_t;

    // this test is for two nodes
    if (nodes != 2)
        return;

    // make a partitioner
    SP<Partitioner> p;

    // database
    Simple_Partitioner::SP_Std_DB db(new Std_DB("partition"));
    UNIT_TEST(db);

    // set data
    {
        db->new_key("num_blocks_i", 1);
        db->new_key("num_blocks_j", 2);
        db->new_key("num_cells_i", 10);
        db->new_key("num_cells_j", 10);
        db->new_key("num_cells_k", 10);
        db->new_key("delta_x", 1.0);
        db->new_key("delta_y", 1.0);
        db->new_key("delta_z", 1.0);
        db->new_key("num_z_blocks", 1);
        db->new_key("num_groups", 1);
    }

    // make the simple partitioner
    {
        // simple partitioner specialization
        p = new Simple_Partitioner(db);
        UNIT_TEST(p);

        if (!p) ITFAILS;
    }
    
    // build the mesh
    p->build();

    // pretend we are on a CFD application node
    if (node == 0) 
    {	
        // make a point map
        Point_Map_t map(p);
	
        // make a field of points
        std::vector<DataType> field_pts = make_field_pts<DataType>();
	
        // add the point field to the map
        for (int i = 0; i < 5; ++i) 
        {
            map.add_point( Point_t(field_pts[3*i], field_pts[3*i+1], 
                                   field_pts[3*i+2], i) );
        }

        // make sure all of the field points got into the map
        if (map.num_nodes() != 5) ITFAILS;
	
        // complete the map
        map.complete();
	
        // check the status of the map
        if (!map.status()) ITFAILS;
	
        // insert some fake data into the map
        {
            // point 1
            map.insert_data(0, "power", 1.35);
	    
            // point 2
            map.insert_data(1, "power", 8.964);
	    
            // point 3
            map.insert_data(2, "power", 0.201);
	    
            // point 4
            map.insert_data(3, "power", 7.775);
	    
            // point 5
            map.insert_data(4, "power", 982.2);
        }
	
        // verify the contents of the map
        {
            // number of nodes in the map
            if (map.num_nodes() != 5)           ITFAILS;
	    
            // unique set of partition ids
            if (map.partitions().size() != 3)   ITFAILS;
            if (map.partitions()[0] != -1)      ITFAILS;
            if (map.partitions()[1] != 0)       ITFAILS;
            if (map.partitions()[2] != 1)       ITFAILS;

            // point 1
            if (map.get_node(0).partition() != 0)                   ITFAILS;
            if ( !soft_equiv(map.get_data(0, "power"), 1.35) )      ITFAILS;
            if ( !soft_equiv(map.get_node(0).point().x(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(0).point().y(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(0).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(0).point().handle() != 0)              ITFAILS;
	    
            // point 2
            if (map.get_node(1).partition() != 0)                   ITFAILS;
            if ( !soft_equiv(map.get_data(1, "power"), 8.964) )     ITFAILS;
            if ( !soft_equiv(map.get_node(1).point().x(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(1).point().y(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(1).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(1).point().handle() != 1)              ITFAILS;
	    
            // point 3
            if (map.get_node(2).partition() != 1)                   ITFAILS;
            if ( !soft_equiv(map.get_data(2, "power"), 0.201) )     ITFAILS;
            if ( !soft_equiv(map.get_node(2).point().x(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(2).point().y(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(2).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(2).point().handle() != 2)              ITFAILS;
	    
            // point 4
            if (map.get_node(3).partition() != 1)                   ITFAILS;
            if ( !soft_equiv(map.get_data(3, "power"), 7.775) )     ITFAILS;
            if ( !soft_equiv(map.get_node(3).point().x(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(3).point().y(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(3).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(3).point().handle() != 3)              ITFAILS;
	    
            // point 5
            if (map.get_node(4).partition() != -1)                  ITFAILS;
            if ( !soft_equiv(map.get_data(4, "power"), 982.2) )     ITFAILS;
            if ( !soft_equiv(map.get_node(4).point().x(), -2.5) )   ITFAILS;
            if ( !soft_equiv(map.get_node(4).point().y(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(4).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(4).point().handle() != 4)              ITFAILS;
        }
    }

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "Two processor Point_Map test ok on " << node;
        ut.passes( m.str() );
    }
}

//---------------------------------------------------------------------------//
// four processor test - 2 I blocks, 2 J blocks
void four_pe_map_test(Parallel_Unit_Test &ut)
{
    typedef int                                         OrdinateType;
    typedef double                                      DataType;
    typedef coupler::Point_Map<OrdinateType, DataType>  Point_Map_t;
    typedef typename Point_Map_t::Point_t               Point_t;

    // this test is for four nodes
    if (nodes != 4)
        return;

    // make a partitioner
    SP<Partitioner> p;

    // database
    Simple_Partitioner::SP_Std_DB db(new Std_DB("partition"));
    UNIT_TEST(db);

    // set data
    {
        db->new_key("num_blocks_i", 2);
        db->new_key("num_blocks_j", 2);
        db->new_key("num_cells_i", 10);
        db->new_key("num_cells_j", 10);
        db->new_key("num_cells_k", 10);
        db->new_key("delta_x", 1.0);
        db->new_key("delta_y", 1.0);
        db->new_key("delta_z", 1.0);
        db->new_key("num_z_blocks", 1);
        db->new_key("num_groups", 1);
    }

    // make the simple partitioner
    {
        // simple partitioner specialization
        p = new Simple_Partitioner(db);
        UNIT_TEST(p);

        if (!p) ITFAILS;
    }
    
    // build the mesh
    p->build();

    // pretend we are on a CFD application node
    if (node == 0) 
    {
        // make a point map
        Point_Map_t map(p);
	
        // make a field of points
        std::vector<DataType> field_pts = make_field_pts<DataType>();
	
        // add the point field to the map
        for (OrdinateType i = 0; i < 5; ++i) 
        {
            map.add_point( Point_t(field_pts[3*i], field_pts[3*i+1], 
                                   field_pts[3*i+2], i) );
        }

        // make sure all of the field points got into the map
        if (map.num_nodes() != 5) ITFAILS;
	
        // complete the map
        map.complete();
	
        // check the status of the map
        if ( !map.status() ) ITFAILS;
	
        // insert some fake data into the map
        {
            // point 1
            map.insert_data(0, "power", 1.35);
	    
            // point 2
            map.insert_data(1, "power", 8.964);
	    
            // point 3
            map.insert_data(2, "power", 0.201);
	    
            // point 4
            map.insert_data(3, "power", 7.775);
	    
            // point 5
            map.insert_data(4, "power", 982.2);
        }

        // verify the contents of the map
        {
            // number of nodes in the map
            if (map.num_nodes() != 5)                   ITFAILS;
	    
            // unique set of partition ids
            if (map.partitions().size() != 5)           ITFAILS;
            if (map.partitions()[0] != -1)              ITFAILS;
            if (map.partitions()[1] != 0)               ITFAILS;
            if (map.partitions()[2] != 1)               ITFAILS;
            if (map.partitions()[3] != 2)               ITFAILS;
            if (map.partitions()[4] != 3)               ITFAILS;
	    
            // point 1
            if (map.get_node(0).partition() != 0)                   ITFAILS;
            if ( !soft_equiv(map.get_data(0, "power"), 1.35) )      ITFAILS;
            if ( !soft_equiv(map.get_node(0).point().x(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(0).point().y(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(0).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(0).point().handle() != 0)              ITFAILS;
	    
            // point 2
            if (map.get_node(1).partition() != 1)                   ITFAILS;
            if ( !soft_equiv(map.get_data(1, "power"), 8.964) )     ITFAILS;

            if ( !soft_equiv(map.get_node(1).point().x(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(1).point().y(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(1).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(1).point().handle() != 1)              ITFAILS;
	    
            // point 3
            if (map.get_node(2).partition() != 2)                   ITFAILS;
            if ( !soft_equiv(map.get_data(2, "power"), 0.201) )     ITFAILS;
            if ( !soft_equiv(map.get_node(2).point().x(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(2).point().y(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(2).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(2).point().handle() != 2)              ITFAILS;
	    
            // point 4
            if (map.get_node(3).partition() != 3)                   ITFAILS;
            if ( !soft_equiv(map.get_data(3, "power"), 7.775) )     ITFAILS;

            if ( !soft_equiv(map.get_node(3).point().x(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(3).point().y(), 7.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(3).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(3).point().handle() != 3)              ITFAILS;
	    
            // point 5
            if (map.get_node(4).partition() != -1)                  ITFAILS;
            if ( !soft_equiv(map.get_data(4, "power"), 982.2) )     ITFAILS;
            if ( !soft_equiv(map.get_node(4).point().x(), -2.5) )   ITFAILS;
            if ( !soft_equiv(map.get_node(4).point().y(), 2.5) )    ITFAILS;
            if ( !soft_equiv(map.get_node(4).point().z(), 2.5) )    ITFAILS;
            if (map.get_node(4).point().handle() != 4)              ITFAILS;
        }
    }

    if (ut.numFails == 0)
    {
        std::ostringstream m;
        m << "Four processor Point_Map test ok on " << node;
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
            one_pe_map_test(ut);
            gpass += ut.numPasses;
            gfail += ut.numFails;
            ut.reset();
        }

        else if (nodes == 2)
        {
            two_pe_map_test(ut);
            gpass += ut.numPasses;
            gfail += ut.numFails;
            ut.reset();

            two_pe_switch_map_test(ut);
            gpass += ut.numPasses;
            gfail += ut.numFails;
            ut.reset();
        }

        else if (nodes == 4)
        {
            four_pe_map_test(ut);
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
        std::cout << "ERROR: While testing tstPoint_Map, " 
                  << err.what()
                  << std::endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstPoint_Map, " 
                  << "An unknown exception was thrown."
                  << std::endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstPoint_Map.cc
//---------------------------------------------------------------------------//
