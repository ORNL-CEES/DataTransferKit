//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/test/tstCoupler_Map_Gen.cc
 * \author Stuart R. Slattery
 * \date   Fri Jun 17 14:55:53 2011
 * \brief  Explicit test for the general coupler map building method.
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
#include "utils/Packing_Utils.hh"
#include "utils/SP.hh"
#include "utils/Vector_Lite.hh"
#include "comm/global.hh"
#include "comm/Parallel_Unit_Test.hh"
#include "release/Release.hh"
#include "mesh_type/Point.hh"
#include "kba_mesh/Partitioner.hh"
#include "kba_mesh/Simple_Partitioner.hh"
#include "database/Std_DB.hh"
#include "../Map_Node.hh"
#include "../Point_Map.hh"
#include "../LG_Indexer.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

using database::Std_DB;

using denovo::Point;
using denovo::SP;

using kba::Partitioner;
using kba::Simple_Partitioner;

using coupler::Point_Map;
using coupler::Map_Node;
using coupler::LG_Indexer;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// make a vector of point objects with globally unique handles
vector<Point<double> > make_points(int node)
{
    vector<Point<double> > points;

    // point 1
    Point<double> p1( 0.0043, 1.322, 5.3365, 4*node+0);
    points.push_back(p1);

    // point 2
    Point<double> p2( 9.56, 7.56, 2.36, 4*node+1);
    points.push_back(p2);

    // point 3
    Point<double> p3( 2.654, 8.225, 5.5556, 4*node+2);
    points.push_back(p3);

    // point 4
    Point<double> p4( 6.5458, 3.748, 7.1645, 4*node+3);
    points.push_back(p4);

    return points;
}

//---------------------------------------------------------------------------//
// simple neutronics class
class Neutronics
{
  public:
    Neutronics() { /* ... */ }
    ~Neutronics() { /* ... */ }
};

//---------------------------------------------------------------------------//
// simple cfd class
class CFD
{
  public:
    CFD() { /* ... */ }
    ~CFD() { /* ... */ }
};

//---------------------------------------------------------------------------//
// make a simple partitioner
SP<Simple_Partitioner> make_partitioner(int num_i_blocks, int num_j_blocks)
{
    // make a partitioner
    SP<Partitioner> p;

    // database
    Simple_Partitioner::SP_Std_DB db(new Std_DB("partition"));

    // set data
    {
        db->new_key("num_blocks_i", num_i_blocks);
        db->new_key("num_blocks_j", num_j_blocks);
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
    }
    
    // build the mesh
    p->build();

    return p;
}

//---------------------------------------------------------------------------//
// point comparator
static bool handle_compare(Point<double> p1, Point<double> p2)
{
    return p1.get_handle() < p2.get_handle();
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
// Explicitly check the code in Coupler::build_map(). We need to do this here
// so we have direct access to the Point_Map object which is not exposed by
// the Coupler in order to verify its contents.
void build_map_test(Parallel_Unit_Test &ut)
{
    typedef vector<int>                  Vec_Int;
    typedef Point<double>                Point_Dbl;
    typedef vector<Point_Dbl>            Vec_Point;
    typedef SP<LG_Indexer>               SP_LG_Indexer;
    typedef SP<Simple_Partitioner>       SP_Simple_Partitioner;
    typedef nemesis::Communicator_t      Communicator_t;
    typedef vector<char>                 Buffer;


    // partitioner setup
    int num_i_blocks, num_j_blocks;
    if (nodes == 1)
    {
	num_i_blocks = 1;
	num_j_blocks = 1;
    }
    else if (nodes == 2)
    {
	num_i_blocks = 2;
	num_j_blocks = 1;
    }
    else if (nodes == 4)
    {
	num_i_blocks = 2;
	num_j_blocks = 2;
    }

    // setup neutronics
    nemesis::set_internal_comm( get_comm_world() );
    Communicator_t comm_neutronics;
    nemesis::split(0, nemesis::nodes()-nemesis::node()-1, comm_neutronics);
    nemesis::set_internal_comm(comm_neutronics);
    Point_Map neutronics_map;
    Vec_Point neutronics_points;
    SP_Simple_Partitioner partitioner = 
	make_partitioner(num_i_blocks, num_j_blocks);
    denovo::SP<Neutronics> neutronics( new Neutronics() );

    // setup cfd
    nemesis::set_internal_comm( get_comm_world() );
    Communicator_t comm_cfd;
    nemesis::split(0, nemesis::node(), comm_cfd);
    nemesis::set_internal_comm(comm_cfd);
    Point_Map cfd_map;
    Vec_Point cfd_points = make_points(nemesis::node());
    denovo::SP<CFD> cfd( new CFD() );

    // make indexers
    nemesis::set_internal_comm( get_comm_world() );
    SP_LG_Indexer indexer_neutronics( 
        new LG_Indexer(get_comm_world(), comm_neutronics, neutronics) );
    SP_LG_Indexer indexer_cfd( 
        new LG_Indexer(get_comm_world(), comm_cfd, cfd) );

    // build the local map

    // set local process id iteration bounds
    int neutronics_begin = 0;
    int neutronics_end = indexer_neutronics->size();
    int external_begin = 0;
    int external_end = indexer_cfd->size();
			       
    // send points from CFD to neutronics
    // create a packer to send points to neutronics
    denovo::Packer p;

    // compute the size of the buffer
    p.compute_buffer_size_mode();
    for (Vec_Point::const_iterator iter = cfd_points.begin(),
			       iter_end = cfd_points.end();
	 iter != iter_end; ++iter)
    {
	p << *iter;
    }
    int buffer_size = p.size();
    Check(buffer_size >= 0);

    // pack the points into the buffer
    Buffer buffer(buffer_size);
    if (buffer_size > 0)
    {
	p.set_buffer(buffer_size, &buffer[0]);
	for (Vec_Point::const_iterator iter = cfd_points.begin(),
				   iter_end = cfd_points.end();
	     iter != iter_end; ++iter)
	{
	    p << *iter;
	}
    }

    // non-blocking send the local points to all neutronics processes
    for (int i = neutronics_begin; i < neutronics_end; ++i)
    {
	// get the global index for neutronics that we are sending to
	int destination = indexer_neutronics->l2g(i);

	// send a message so neutronics knows what is coming
	nemesis::send_async(&buffer_size, 1, destination);

	// send the actual buffer
	if (buffer_size > 0)
	{
	    nemesis::send_async(&buffer[0], buffer_size, destination);
	}
    }

    // neutronics gets points from cfd and makes its map
    // every external process sends neutronics it's points
    for (int i = external_begin; i < external_end; ++i)
    {
	// get the global index for external that we are receiving from
	int source = indexer_cfd->l2g(i);

	// receive number of points from the external application
	int buffer_size;
	nemesis::receive(&buffer_size, 1, source);

	// if the buffer is not empty
	if (buffer_size > 0)
	{
	    // receive the buffer from the external application
	    Buffer buffer(buffer_size);
	    nemesis::receive(&buffer[0], buffer_size, source);

	    // compute the number of point objects in the buffer
	    int num_points = buffer_size / sizeof(Point_Dbl);

	    // unpack the buffer
	    denovo::Unpacker u;
	    u.set_buffer(buffer_size, &buffer[0]);
	    for (int j = 0; j < num_points; ++j)
	    {
		Point_Dbl point(0.0, 0.0, 0.0, 0);
		u >> point;

		// see if the point is in this process' spatial domain
		denovo::Vector_Lite<int, 3> dummy_vector;
		if ( partitioner->get_mesh()->find_cell(
			 point.get_coords(), dummy_vector) )
		{
		    // Make a new map node
		    Map_Node map_node(point.get_handle(), source, point);

		    // add it to the map
		    neutronics_map.add_map_node(map_node);
		}
	    }
	}
    }
	    
    // complete the neutronics mapping
    neutronics_map.complete();

    // barrier after the neutronics map generation
    nemesis::global_barrier();

    // neutronics sends points back to cfd
    // get the unique partition id vector for this neutronics node
    Vec_Int partitions = neutronics_map.get_partitions();

    // get the number of points this neutronics node is sending to each
    // external process
    Vec_Int num_send_points = neutronics_map.get_partition_pts();

    // loop through the map nodes for each destination process and send
    // the points back in a buffer
    int node_counter = 0;
    for (int i = 0; i < partitions.size(); ++i)
    {
	// create a Packer
	denovo::Packer p;

	// compute the size of the buffer
	p.compute_buffer_size_mode();
	for (int j = 0; j < num_send_points[i]; ++j)
	{
	    p << 
		neutronics_map.get_all_nodes()[node_counter].get_node_point();
	    ++node_counter;
	}

	// reset the node counter
	node_counter -= num_send_points[i];

	// create and pack the buffer
	int buffer_size = p.size();
	Buffer buffer(buffer_size);
	p.set_buffer(buffer_size, &buffer[0]);
	for (int j = 0; j < num_send_points[i]; ++j)
	{
	    p << neutronics_map.get_all_nodes()[node_counter].get_node_point();
	    ++node_counter;
	}

	// get the global index to external that we are sending to
	int destination = partitions[i];

	// send the size of the buffer
	nemesis::send_async(&buffer_size, 1, destination);

	// send the buffer
	nemesis::send_async(&buffer[0], buffer_size, destination);
    }

    // cfd gets the points and finishes its map
    // get a buffer from all neutronics processes
    for (int i = neutronics_begin; i < neutronics_end; ++i)
    {
	// get the global index to neutronics that we are receiving from
	int source = indexer_neutronics->l2g(i);
	    
	// receive the buffer size
	int buffer_size;
	nemesis::receive(&buffer_size, 1, source);

	if (buffer_size > 0)
	{
	    // create a buffer
	    Buffer buffer(buffer_size);

	    // receive the message
	    nemesis::receive(&buffer[0], buffer_size, source);

	    // compute the number of points in the buffer
	    int num_points = buffer_size / sizeof(Point_Dbl);

	    // unpack the buffer
	    denovo::Unpacker u;
	    u.set_buffer(buffer_size, &buffer[0]);
	    for (int j = 0; j < num_points; ++j)
	    {
		// get the point
		Point_Dbl dummy(0.0, 0.0, 0.0, 0);
		u >> dummy;

		// find the local point with the matching handle
		Vec_Point::iterator it = std::lower_bound(cfd_points.begin(), 
							  cfd_points.end(),
							  dummy,
							  &handle_compare);

		// add a node to the map for the point
		Map_Node map_node(it->get_handle(), source, *it);
		cfd_map.add_map_node(map_node);
	    }
	}
    }

    // complete the external application mapping
    cfd_map.complete();

    // barrier after the external map generation
    nemesis::global_barrier();

    // check both maps
    if (nodes == 1)
    {
	UNIT_TEST(neutronics_map.get_num_nodes() == 4);
	UNIT_TEST(neutronics_map.get_partitions().size() == 1);
	UNIT_TEST(neutronics_map.get_partitions()[0] == 0);
	UNIT_TEST(neutronics_map.get_partition_pts().size() == 1);
	UNIT_TEST(neutronics_map.get_partition_pts()[0] == 4);
	UNIT_TEST(neutronics_map.get_map_node(0).get_node_point().x() 
		  == 0.0043);
	UNIT_TEST(neutronics_map.get_map_node(1).get_node_point().x() 
		  == 9.56);
	UNIT_TEST(neutronics_map.get_map_node(2).get_node_point().x() 
		  == 2.654);
	UNIT_TEST(neutronics_map.get_map_node(3).get_node_point().x() 
		  == 6.5458);
	UNIT_TEST(neutronics_map.get_all_nodes()[0].get_node_point().x() 
		  == 0.0043);
	UNIT_TEST(neutronics_map.get_all_nodes()[1].get_node_point().x() 
		  == 9.56);
	UNIT_TEST(neutronics_map.get_all_nodes()[2].get_node_point().x() 
		  == 2.654);
	UNIT_TEST(neutronics_map.get_all_nodes()[3].get_node_point().x() 
		  == 6.5458);

	UNIT_TEST(cfd_map.get_num_nodes() == 4);
	UNIT_TEST(cfd_map.get_partitions().size() == 1);
	UNIT_TEST(cfd_map.get_partitions()[0] == 0);
	UNIT_TEST(cfd_map.get_partition_pts().size() == 1);
	UNIT_TEST(cfd_map.get_partition_pts()[0] == 4);
	UNIT_TEST(cfd_map.get_map_node(0).get_node_point().x() 
		  == 0.0043);
	UNIT_TEST(cfd_map.get_map_node(1).get_node_point().x() 
		  == 9.56);
	UNIT_TEST(cfd_map.get_map_node(2).get_node_point().x() 
		  == 2.654);
	UNIT_TEST(cfd_map.get_map_node(3).get_node_point().x() 
		  == 6.5458);
	UNIT_TEST(cfd_map.get_all_nodes()[0].get_node_point().x() 
		  == 0.0043);
	UNIT_TEST(cfd_map.get_all_nodes()[1].get_node_point().x() 
		  == 9.56);
	UNIT_TEST(cfd_map.get_all_nodes()[2].get_node_point().x() 
		  == 2.654);
	UNIT_TEST(cfd_map.get_all_nodes()[3].get_node_point().x() 
		  == 6.5458);
    }
    if (nodes == 2)
    {
	if (node == 0)
	{
	    UNIT_TEST(neutronics_map.get_num_nodes() == 4);
	    UNIT_TEST(neutronics_map.get_partitions().size() == 2);
	    UNIT_TEST(neutronics_map.get_partitions()[0] == 0);
	    UNIT_TEST(neutronics_map.get_partitions()[1] == 1);	    
	    UNIT_TEST(neutronics_map.get_partition_pts().size() == 2);
	    UNIT_TEST(neutronics_map.get_partition_pts()[0] == 2);
	    UNIT_TEST(neutronics_map.get_partition_pts()[1] == 2);
	    UNIT_TEST(neutronics_map.get_map_node(1).get_node_point().x()
		      == 9.56);
	    UNIT_TEST(neutronics_map.get_map_node(3).get_node_point().x() 
		      == 6.5458);
	    UNIT_TEST(neutronics_map.get_map_node(5).get_node_point().x()
		      == 9.56);
	    UNIT_TEST(neutronics_map.get_map_node(7).get_node_point().x() 
		      == 6.5458);

	    UNIT_TEST(cfd_map.get_num_nodes() == 4);
	    UNIT_TEST(cfd_map.get_partitions().size() == 2);
	    UNIT_TEST(cfd_map.get_partitions()[0] == 0);
	    UNIT_TEST(cfd_map.get_partitions()[1] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts().size() == 2);
	    UNIT_TEST(cfd_map.get_partition_pts()[0] == 2);
	    UNIT_TEST(cfd_map.get_partition_pts()[1] == 2);
	    UNIT_TEST(cfd_map.get_map_node(0).get_node_point().x() 
		      == 0.0043);
	    UNIT_TEST(cfd_map.get_map_node(1).get_node_point().x() 
		      == 9.56);
	    UNIT_TEST(cfd_map.get_map_node(2).get_node_point().x() 
		      == 2.654);
	    UNIT_TEST(cfd_map.get_map_node(3).get_node_point().x() 
		      == 6.5458);
	}
	if (node == 1)
	{
	    UNIT_TEST(neutronics_map.get_num_nodes() == 4);
	    UNIT_TEST(neutronics_map.get_partitions().size() == 2);
	    UNIT_TEST(neutronics_map.get_partitions()[0] == 0);
	    UNIT_TEST(neutronics_map.get_partitions()[1] == 1);	    
	    UNIT_TEST(neutronics_map.get_partition_pts().size() == 2);
	    UNIT_TEST(neutronics_map.get_partition_pts()[0] == 2);
	    UNIT_TEST(neutronics_map.get_partition_pts()[1] == 2);
	    UNIT_TEST(neutronics_map.get_map_node(0).get_node_point().x() 
		      == 0.0043);
	    UNIT_TEST(neutronics_map.get_map_node(2).get_node_point().x() 
		      == 2.654);
	    UNIT_TEST(neutronics_map.get_map_node(4).get_node_point().x() 
		      == 0.0043);
	    UNIT_TEST(neutronics_map.get_map_node(6).get_node_point().x() 
		      == 2.654);

	    UNIT_TEST(cfd_map.get_num_nodes() == 4);
	    UNIT_TEST(cfd_map.get_partitions().size() == 2);
	    UNIT_TEST(cfd_map.get_partitions()[0] == 0);
	    UNIT_TEST(cfd_map.get_partitions()[1] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts().size() == 2);
	    UNIT_TEST(cfd_map.get_partition_pts()[0] == 2);
	    UNIT_TEST(cfd_map.get_partition_pts()[1] == 2);
	    UNIT_TEST(cfd_map.get_map_node(4).get_node_point().x() 
		      == 0.0043);
	    UNIT_TEST(cfd_map.get_map_node(5).get_node_point().x() 
		      == 9.56);
	    UNIT_TEST(cfd_map.get_map_node(6).get_node_point().x() 
		      == 2.654);
	    UNIT_TEST(cfd_map.get_map_node(7).get_node_point().x() 
		      == 6.5458);
	}
    }
    if (nodes == 4)
    {
	if (node == 0)
	{
	    UNIT_TEST(neutronics_map.get_num_nodes() == 4);
	    UNIT_TEST(neutronics_map.get_partitions().size() == 4);
	    UNIT_TEST(neutronics_map.get_partitions()[0] == 0);
	    UNIT_TEST(neutronics_map.get_partitions()[1] == 1);
	    UNIT_TEST(neutronics_map.get_partitions()[2] == 2);
	    UNIT_TEST(neutronics_map.get_partitions()[3] == 3);
	    UNIT_TEST(neutronics_map.get_partition_pts().size() == 4);
	    UNIT_TEST(neutronics_map.get_partition_pts()[0] == 1);
	    UNIT_TEST(neutronics_map.get_partition_pts()[1] == 1);
	    UNIT_TEST(neutronics_map.get_partition_pts()[2] == 1);
	    UNIT_TEST(neutronics_map.get_partition_pts()[3] == 1);
	    UNIT_TEST(neutronics_map.get_map_node(1).get_node_point().x() 
		      == 9.56);
	    UNIT_TEST(neutronics_map.get_map_node(5).get_node_point().x() 
		      == 9.56);
	    UNIT_TEST(neutronics_map.get_map_node(9).get_node_point().x() 
		      == 9.56);
	    UNIT_TEST(neutronics_map.get_map_node(13).get_node_point().x() 
		      == 9.56);

	    UNIT_TEST(cfd_map.get_num_nodes() == 4);
	    UNIT_TEST(cfd_map.get_partitions().size() == 4);
	    UNIT_TEST(cfd_map.get_partitions()[0] == 0);
	    UNIT_TEST(cfd_map.get_partitions()[1] == 1);
	    UNIT_TEST(cfd_map.get_partitions()[2] == 2);
	    UNIT_TEST(cfd_map.get_partitions()[3] == 3);
	    UNIT_TEST(cfd_map.get_partition_pts().size() == 4);
	    UNIT_TEST(cfd_map.get_partition_pts()[0] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts()[1] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts()[2] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts()[3] == 1);
	    UNIT_TEST(cfd_map.get_map_node(0).get_node_point().x() 
		      == 0.0043);
	    UNIT_TEST(cfd_map.get_map_node(1).get_node_point().x() 
		      == 9.56);
	    UNIT_TEST(cfd_map.get_map_node(2).get_node_point().x() 
		      == 2.654);
	    UNIT_TEST(cfd_map.get_map_node(3).get_node_point().x() 
		      == 6.5458);
	}
	if (node == 1)
	{
	    UNIT_TEST(neutronics_map.get_num_nodes() == 4);
	    UNIT_TEST(neutronics_map.get_partitions().size() == 4);
	    UNIT_TEST(neutronics_map.get_partitions()[0] == 0);
	    UNIT_TEST(neutronics_map.get_partitions()[1] == 1);
	    UNIT_TEST(neutronics_map.get_partitions()[2] == 2);
	    UNIT_TEST(neutronics_map.get_partitions()[3] == 3);
	    UNIT_TEST(neutronics_map.get_partition_pts().size() == 4);
	    UNIT_TEST(neutronics_map.get_partition_pts()[0] == 1);
	    UNIT_TEST(neutronics_map.get_partition_pts()[1] == 1);
	    UNIT_TEST(neutronics_map.get_partition_pts()[2] == 1);
	    UNIT_TEST(neutronics_map.get_partition_pts()[3] == 1);
	    UNIT_TEST(neutronics_map.get_map_node(2).get_node_point().x() 
		      == 2.654);
	    UNIT_TEST(neutronics_map.get_map_node(6).get_node_point().x() 
		      == 2.654);
	    UNIT_TEST(neutronics_map.get_map_node(10).get_node_point().x() 
		      == 2.654);
	    UNIT_TEST(neutronics_map.get_map_node(14).get_node_point().x() 
		      == 2.654);

	    UNIT_TEST(cfd_map.get_num_nodes() == 4);
	    UNIT_TEST(cfd_map.get_partitions().size() == 4);
	    UNIT_TEST(cfd_map.get_partitions()[0] == 0);
	    UNIT_TEST(cfd_map.get_partitions()[1] == 1);
	    UNIT_TEST(cfd_map.get_partitions()[2] == 2);
	    UNIT_TEST(cfd_map.get_partitions()[3] == 3);
	    UNIT_TEST(cfd_map.get_partition_pts().size() == 4);
	    UNIT_TEST(cfd_map.get_partition_pts()[0] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts()[1] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts()[2] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts()[3] == 1);
	    UNIT_TEST(cfd_map.get_map_node(4).get_node_point().x() 
		      == 0.0043);
	    UNIT_TEST(cfd_map.get_map_node(5).get_node_point().x() 
		      == 9.56);
	    UNIT_TEST(cfd_map.get_map_node(6).get_node_point().x() 
		      == 2.654);
	    UNIT_TEST(cfd_map.get_map_node(7).get_node_point().x() 
		      == 6.5458);
	}
	if (node == 2)
	{
	    UNIT_TEST(neutronics_map.get_num_nodes() == 4);
	    UNIT_TEST(neutronics_map.get_partitions().size() == 4);
	    UNIT_TEST(neutronics_map.get_partitions()[0] == 0);
	    UNIT_TEST(neutronics_map.get_partitions()[1] == 1);
	    UNIT_TEST(neutronics_map.get_partitions()[2] == 2);
	    UNIT_TEST(neutronics_map.get_partitions()[3] == 3);
	    UNIT_TEST(neutronics_map.get_partition_pts().size() == 4);
	    UNIT_TEST(neutronics_map.get_partition_pts()[0] == 1);
	    UNIT_TEST(neutronics_map.get_partition_pts()[1] == 1);
	    UNIT_TEST(neutronics_map.get_partition_pts()[2] == 1);
	    UNIT_TEST(neutronics_map.get_partition_pts()[3] == 1);
	    UNIT_TEST(neutronics_map.get_map_node(3).get_node_point().x()
		      == 6.5458);
	    UNIT_TEST(neutronics_map.get_map_node(7).get_node_point().x()
		      == 6.5458);
	    UNIT_TEST(neutronics_map.get_map_node(11).get_node_point().x() 
		      == 6.5458);
	    UNIT_TEST(neutronics_map.get_map_node(15).get_node_point().x() 
		      == 6.5458);

	    UNIT_TEST(cfd_map.get_num_nodes() == 4);
	    UNIT_TEST(cfd_map.get_partitions().size() == 4);
	    UNIT_TEST(cfd_map.get_partitions()[0] == 0);
	    UNIT_TEST(cfd_map.get_partitions()[1] == 1);
	    UNIT_TEST(cfd_map.get_partitions()[2] == 2);
	    UNIT_TEST(cfd_map.get_partitions()[3] == 3);
	    UNIT_TEST(cfd_map.get_partition_pts().size() == 4);
	    UNIT_TEST(cfd_map.get_partition_pts()[0] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts()[1] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts()[2] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts()[3] == 1);
	    UNIT_TEST(cfd_map.get_map_node(8).get_node_point().x() 
		      == 0.0043);
	    UNIT_TEST(cfd_map.get_map_node(9).get_node_point().x() 
		      == 9.56);
	    UNIT_TEST(cfd_map.get_map_node(10).get_node_point().x() 
		      == 2.654);
	    UNIT_TEST(cfd_map.get_map_node(11).get_node_point().x() 
		      == 6.5458);
	}
	if (node == 3)
	{
	    UNIT_TEST(neutronics_map.get_num_nodes() == 4);
	    UNIT_TEST(neutronics_map.get_partitions().size() == 4);
	    UNIT_TEST(neutronics_map.get_partitions()[0] == 0);
	    UNIT_TEST(neutronics_map.get_partitions()[1] == 1);
	    UNIT_TEST(neutronics_map.get_partitions()[2] == 2);
	    UNIT_TEST(neutronics_map.get_partitions()[3] == 3);
	    UNIT_TEST(neutronics_map.get_partition_pts().size() == 4);
	    UNIT_TEST(neutronics_map.get_partition_pts()[0] == 1);
	    UNIT_TEST(neutronics_map.get_partition_pts()[1] == 1);
	    UNIT_TEST(neutronics_map.get_partition_pts()[2] == 1);
	    UNIT_TEST(neutronics_map.get_partition_pts()[3] == 1);
	    UNIT_TEST(neutronics_map.get_map_node(0).get_node_point().x() 
		      == 0.0043);
	    UNIT_TEST(neutronics_map.get_map_node(4).get_node_point().x() 
		      == 0.0043);
	    UNIT_TEST(neutronics_map.get_map_node(8).get_node_point().x() 
		      == 0.0043);
	    UNIT_TEST(neutronics_map.get_map_node(12).get_node_point().x()
		      == 0.0043);

	    UNIT_TEST(cfd_map.get_num_nodes() == 4);
	    UNIT_TEST(cfd_map.get_partitions().size() == 4);
	    UNIT_TEST(cfd_map.get_partitions()[0] == 0);
	    UNIT_TEST(cfd_map.get_partitions()[1] == 1);
	    UNIT_TEST(cfd_map.get_partitions()[2] == 2);
	    UNIT_TEST(cfd_map.get_partitions()[3] == 3);
	    UNIT_TEST(cfd_map.get_partition_pts().size() == 4);
	    UNIT_TEST(cfd_map.get_partition_pts()[0] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts()[1] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts()[2] == 1);
	    UNIT_TEST(cfd_map.get_partition_pts()[3] == 1);
	    UNIT_TEST(cfd_map.get_map_node(12).get_node_point().x() 
		      == 0.0043);
	    UNIT_TEST(cfd_map.get_map_node(13).get_node_point().x() 
		      == 9.56);
	    UNIT_TEST(cfd_map.get_map_node(14).get_node_point().x() 
		      == 2.654);
	    UNIT_TEST(cfd_map.get_map_node(15).get_node_point().x() 
		      == 6.5458);
	}
    }

    if (ut.numFails == 0)
    {
	ostringstream m;
	m << "Coupler map generation test OK on " << node;
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
            build_map_test(ut);
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
        std::cout << "ERROR: While testing tstCoupler_Map_Gen, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstCoupler_Map_Gen, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstCoupler_Map_Gen.cc
//---------------------------------------------------------------------------//
