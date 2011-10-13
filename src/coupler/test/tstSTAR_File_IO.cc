//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/test/tstSTAR_Parser.cc
 * \author Gregory Davidson
 * \date   Tue Jun 07 16:19:28 2011
 * \brief  Tests the STAR_Parser class.
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
#include "utils/SP.hh"
#include "database/Std_DB.hh"
#include "release/Release.hh"
#include "material/Mat_DB.hh"
#include "kba_mesh/Partitioner.hh"
#include "kba_mesh/Simple_Partitioner.hh"
#include "kba_equations/SC_Equations.hh"
#include "kba/Erg_Set_Communicator.hh"
#include "../Point_Map.hh"
#include "../STAR_Data.hh"
#include "../STAR_File_IO.hh"

using namespace std;
using nemesis::Parallel_Unit_Test;
using nemesis::soft_equiv;

using denovo::SP;

using database::Std_DB;
using kba::Simple_Partitioner;
using kba::Partitioner;

using coupler::Point_Map;

int node  = 0;
int nodes = 0;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// make STAR points
template<class OrdinateType_T, class DataType_T>
vector< denovo::Point<OrdinateType_T, DataType_T> > make_field_pts()
{
    typedef OrdinateType_T                          OrdinateType;
    typedef DataType_T                              DataType;
    typedef denovo::Point<OrdinateType, DataType>   Point_t;
    typedef std::vector<Point_t>                    Vec_Point;

    Vec_Point pts;

    // first point
    Point_t p1(2.5, 2.5, 2.5, 1);
    pts.push_back(p1);

    // second point
    Point_t p2(7.5, 2.5, 2.5, 2);
    pts.push_back(p2);

    // third point
    Point_t p3(2.5, 7.5, 2.5, 3);
    pts.push_back(p3);

    // fourth point
    Point_t p4(7.5, 7.5, 2.5, 4);
    pts.push_back(p4);

    // fifth point (outside of global mesh)
    Point_t p5(1.7, 2.7, 3.7, 5);
    pts.push_back(p5);

    return pts;
}

//---------------------------------------------------------------------------//
// Create the database
denovo::SP<database::Std_DB> create_db(int num_blocks_I, int num_blocks_J, 
                                       int num_sets,
                                       Parallel_Unit_Test &ut)
{
    typedef denovo::SP<database::Std_DB>        SP_Std_DB;

    UNIT_TEST( num_blocks_I * num_blocks_J * num_sets == nemesis::nodes() );

    SP_Std_DB db( new database::Std_DB("test_db") );
    UNIT_TEST(db);
    {
        db->new_key("num_blocks_i", num_blocks_I);
        db->new_key("num_blocks_j", num_blocks_J);
        db->new_key("num_z_blocks", 1);
        db->new_key("num_sets", num_sets);

        db->new_key("num_groups", 4);
        db->new_key("Pn_order", 0);
        db->new_key("downscatter", false);

        db->new_key("num_cells_i", 2);
        db->new_key("num_cells_j", 2);
        db->new_key("num_cells_k", 2);
        
        db->new_key("delta_x", 1.0);
        db->new_key("delta_y", 1.0);
        db->new_key("delta_z", 1.0);
    }

    // assign number of moments
    denovo::Moment_Order::set_pn_order(1);

    return db;
}


//---------------------------------------------------------------------------//
// Create the erg_set_communicator
template<class Equations>
denovo::SP< kba::Erg_Set_Communicator<Equations> >
create_erg_comm(denovo::SP<database::Std_DB> db, 
                denovo::SP<denovo::Mat_DB> mat,
                Parallel_Unit_Test &ut)
{
    typedef kba::Erg_Set_Communicator<Equations>        Erg_Set_Comm;
    typedef denovo::SP<Erg_Set_Comm>                    SP_Erg_Set_Comm;
    typedef denovo::SP<kba::LG_Indexer>                 SP_Indexer;

    Require (db);

    // Build the partitioner
    kba::Simple_Partitioner p(db);
    p.build();
    // Get the indexer
    SP_Indexer indexer = p.get_indexer();
    UNIT_TEST(indexer);

    // Create the Erg_Set_Comm
    SP_Erg_Set_Comm erg_comm( new Erg_Set_Comm(*indexer, false) );
    UNIT_TEST(erg_comm);

    // Decompose
    erg_comm->decompose(mat, false);

    // Return
    return erg_comm;
}

//---------------------------------------------------------------------------//
// Create a mesh
denovo::SP<kba::Mesh_3D> make_mesh(denovo::SP<database::Std_DB> db,
                                   Parallel_Unit_Test &ut)
{
    typedef denovo::SP<kba::Mesh_3D>        SP_Mesh;

    kba::Simple_Partitioner p(db);
    p.build();
    
    SP_Mesh mesh = p.get_mesh();
    UNIT_TEST(mesh);

    return mesh;
}

//---------------------------------------------------------------------------//
// Create a material database
denovo::SP<denovo::Mat_DB> make_mat(denovo::SP<database::Std_DB> db,
                                    denovo::SP<kba::Mesh_3D> mesh,
                                    Parallel_Unit_Test &ut)
{
    typedef denovo::Mat_DB                      Mat_DB;
    typedef denovo::SP<Mat_DB>                  SP_Mat_DB;

    SP_Mat_DB mat( new denovo::Mat_DB(mesh->num_cells()) );
    UNIT_TEST(mat);
    mat->set_num(1, 4);

    Mat_DB::Vec_XS xs(4);

    for(int g = 0; g < 4; ++g)
    {
        xs[g] = new Mat_DB::XS(db, g);
        UNIT_TEST(xs[g]);

        xs[g]->sigma() = 1.0;
    }

    xs[0]->sigma_s(0, 0) = 0.5;

    xs[1]->sigma_s(0, 0) = 0.5; xs[1]->sigma_s(1, 0) = 0.2;

    xs[2]->sigma_s(0, 0) = 0.5; xs[2]->sigma_s(1, 0) = 0.2;
    xs[2]->sigma_s(2, 0) = 0.1; xs[2]->sigma_s(3, 0) = 0.1;

    xs[3]->sigma_s(0, 0) = 0.5; xs[3]->sigma_s(1, 0) = 0.2;
    xs[3]->sigma_s(2, 0) = 0.1; xs[3]->sigma_s(3, 0) = 0.1;

    for(int g = 0; g < 4; ++g)
    {
        xs[g]->complete();
        mat->assign(xs[g], 0);
    }

    mat->assign(0);

    return mat;
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
void test_geometry_read(int num_blocks_I, int num_blocks_J, int num_sets,
                        const std::string &filename,
                        Parallel_Unit_Test &ut)
{
    typedef kba::SC_Equations                       SC_Equations;
    typedef denovo::SP<database::Std_DB>            SP_Std_DB;
    typedef kba::Erg_Set_Communicator<SC_Equations> Erg_Set_Comm;
    typedef denovo::SP<Erg_Set_Comm>                SP_Erg_Set_Comm;
    typedef denovo::SP<kba::Mesh_3D>                SP_Mesh;
    typedef denovo::SP<denovo::Mat_DB>              SP_Mat_DB;

    typedef coupler::STAR_File_IO<int, double>      STAR_File_IO_t;
    typedef typename STAR_File_IO_t::Vec_Point      Vec_Point;

    // Create the database
    SP_Std_DB db = create_db(num_blocks_I, num_blocks_J, num_sets, ut);
    UNIT_TEST(db);

    // Create the mesh
    SP_Mesh mesh = make_mesh(db, ut);
    Ensure(mesh);

    // Create the material database
    SP_Mat_DB mat = make_mat(db, mesh, ut);
    Ensure(mat);

    // Create the Erg_Set_Comm
    SP_Erg_Set_Comm erg_comm = create_erg_comm<SC_Equations>(db, mat, ut);
    UNIT_TEST(erg_comm);

    // Create the STAR_Parser
    STAR_File_IO_t star_file_io( erg_comm->set_const_comm(), 
                                 erg_comm->block_const_comm() );

    // Read in the test geometry file.
    star_file_io.read_geometry_file("test_star_geom.inp");
    Vec_Point points = star_file_io.get_point_data();

    if(nemesis::node() == 0)
    {
        UNIT_TEST(points.size() == 10);

        // We'll just check the first three lines...

        // First line...
        UNIT_TEST( soft_equiv(points[0].x(), -1.522832e-003) );
        UNIT_TEST( soft_equiv(points[0].y(), 1.239951e-002) );
        UNIT_TEST( soft_equiv(points[0].z(), 1.000004e-002) );

        // Second line...
        UNIT_TEST( soft_equiv(points[1].x(), -1.419053e-003) );
        UNIT_TEST( soft_equiv(points[1].y(), 1.201220e-002) );
        UNIT_TEST( soft_equiv(points[1].z(), 1.000000e-002) );

        // Third line...
        UNIT_TEST( soft_equiv(points[2].x(), -1.522832e-003) );
        UNIT_TEST( soft_equiv(points[2].y(), 1.280048e-002) );
        UNIT_TEST( soft_equiv(points[2].z(), 1.000000e-002) );
    }

    if (ut.numFails == 0)
    {
        ut.passes("Reading STAR Geometry File test ok");
    }
}
        
//---------------------------------------------------------------------------//
// Write the STAR power file
void test_power_write(int num_blocks_I, int num_blocks_J, int num_sets, 
                      const std::string &filename, Parallel_Unit_Test &ut)
{
    typedef kba::SC_Equations                       SC_Equations;
    typedef denovo::SP<database::Std_DB>            SP_Std_DB;
    typedef kba::Erg_Set_Communicator<SC_Equations> Erg_Set_Comm;
    typedef denovo::SP<Erg_Set_Comm>                SP_Erg_Set_Comm;
    typedef denovo::SP<kba::Mesh_3D>                SP_Mesh;
    typedef denovo::SP<denovo::Mat_DB>              SP_Mat_DB;

    typedef int                                             OrdinateType;
    typedef double                                          DataType;
    typedef coupler::STAR_File_IO<OrdinateType, DataType>   STAR_File_IO_t;
    typedef typename STAR_File_IO_t::Point_t                Point_t;
    typedef typename STAR_File_IO_t::Vec_Point              Vec_Point;
    typedef typename STAR_File_IO_t::Point_Map_t            Point_Map_t;
    typedef typename STAR_File_IO_t::SP_Point_Map           SP_Point_Map;

    typedef kba::Simple_Partitioner               Simple_Partitioner;
    typedef denovo::SP<Simple_Partitioner>        SP_Simple_Partitioner;

    // Create the database
    SP_Std_DB db = create_db(num_blocks_I, num_blocks_J, num_sets, ut);
    UNIT_TEST(db);

    // Create the partitioner
    SP_Simple_Partitioner partitioner( new Simple_Partitioner(db) );
    Ensure (partitioner);
    partitioner->build();

    // Create the mesh
    SP_Mesh mesh = make_mesh(db, ut);
    Ensure(mesh);

    // Create the material database
    SP_Mat_DB mat = make_mat(db, mesh, ut);
    Ensure(mat);

    // Create the Erg_Set_Comm
    SP_Erg_Set_Comm erg_comm = create_erg_comm<SC_Equations>(db, mat, ut);
    UNIT_TEST(erg_comm);

    // Create the STAR_Parser
    STAR_File_IO_t star_file_io( erg_comm->set_const_comm(), 
                                 erg_comm->block_const_comm() );

    // Create (on each proc) a set of points
    Vec_Point point_vec;
    for(OrdinateType i = 0; i < 10; ++i)
    {
        DataType x = nemesis::node() + i;
        DataType y = 2.0 * (nemesis::node() + i);
        DataType z = 3.0 * (nemesis::node() + i);

        point_vec.push_back( Point_t(x, y, z, i) );
    }

    // Create a map
    SP_Point_Map point_map( new Point_Map_t(partitioner) );
    Ensure (point_map);

    // Put all the points into the map
    for(Vec_Point::const_iterator iter = point_vec.begin(), 
                              iter_end = point_vec.end(); 
        iter != iter_end; ++iter)
    {
        point_map->add_point(*iter);
    }
    
    // Complete the mapping
    point_map->complete();

    // Add powers to the nodes on each proc
    for(Point_Map_t::iterator iter = point_map->begin(), 
                        iter_end = point_map->end();
        iter != iter_end; ++iter)
    {
        // Make a power
        DataType power = 10.0 * ( nemesis::node() + iter->point().handle() );

        // Register the power
        iter->register_data("power", power);
    }

    // Write the power file
    star_file_io.write_power_file(filename, point_map);

    // Now, read the power file back in (on Node 0) and see if it is right.
    if(nemesis::node() == 0)
    {
        std::ifstream input(filename.c_str());

        // Skip the header
        input.ignore(std::numeric_limits<int>::max(), '\n');

        // Loop over procs
        for(OrdinateType p = 0; p < nemesis::node(); ++p)
        {
            // Loop over data
            for(OrdinateType i = 0; i < 10; ++i)
            {
                DataType power, x, y, z;
                input >> power >> x >> y >> z;

                UNIT_TEST( soft_equiv(power, 10.0 * (p+i) ) );
                UNIT_TEST( soft_equiv(x, 1.0 * (p+i) ) );
                UNIT_TEST( soft_equiv(y, 2.0 * (p+i) ) );
                UNIT_TEST( soft_equiv(z, 3.0 * (p+i) ) );
            }
        }
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Writing STAR Power File test ok on node " << nemesis::node();
        ut.passes(m.str());
    }
}
    
//---------------------------------------------------------------------------//
// serial map building test
void serial_map_build_test(Parallel_Unit_Test &ut)
{
    // typedefs
    typedef int                                             OrdinateType;
    typedef double                                          DataType;
    typedef coupler::STAR_File_IO<OrdinateType, DataType>   STAR_File_IO_t;
    typedef typename STAR_File_IO_t::Point_Map_t            Point_Map_t;
    typedef typename Point_Map_t::Point_t                   Point_t; 
    typedef SP<Point_Map_t>                                 SP_Point_Map;
    typedef typename STAR_File_IO_t::Vec_Point              Vec_Point;

    // this test is serial
    if (nodes > 1)
        return;

    // make a partitioner
    SP<Partitioner> p;

    // database
    Simple_Partitioner::SP_Std_DB db(new Std_DB("partition"));
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
    }

    // build the mesh
    p->build();

    // make a STAR_File_IO object
    STAR_File_IO_t star_io( get_comm_world(), get_comm_world() );

    // build Star data on node 0
    if (node == 0)
    {
        // make some fake Star data
        Vec_Point point_vec = make_field_pts<OrdinateType, DataType>();

        // assign it to the file_io object
        star_io.set_point_data(point_vec);
    }

    // build the local maps on all nodes
    SP_Point_Map map = star_io.build_map(p);

    // check the local map data structure
    if (node == 0) 
    {
        // number of nodes in the map
        if ( map->num_nodes() != 5 )        ITFAILS;

        // unique set of partition ids
        if (map->partitions().size() != 1)  ITFAILS;
        if (map->partitions()[0] != 0)      ITFAILS;

        // point 1
        if (map->get_node(0).partition() != 0)  ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().x(), 2.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().y(), 2.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().z(), 2.5) )   ITFAILS;
        if (map->get_node(0).point().handle() != 1)             ITFAILS;
	    
        // point 2
        if (map->get_node(1).partition() != 0)  ITFAILS;

        if ( !soft_equiv(map->get_node(1).point().x(), 7.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(1).point().y(), 2.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(1).point().z(), 2.5) )   ITFAILS;
        if (map->get_node(1).point().handle() != 2)             ITFAILS;
	    
        // point 3
        if (map->get_node(2).partition() != 0)  ITFAILS;

        if ( !soft_equiv(map->get_node(2).point().x(), 2.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(2).point().y(), 7.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(2).point().z(), 2.5) )   ITFAILS;
        if ( map->get_node(2).point().handle() != 3)            ITFAILS;
	    
        // point 4
        if (map->get_node(3).partition() != 0)  ITFAILS;

        if ( !soft_equiv(map->get_node(3).point().x(), 7.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(3).point().y(), 7.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(3).point().z(), 2.5) )   ITFAILS;
        if (map->get_node(3).point().handle() != 4)             ITFAILS;

        // point 5
        if (map->get_node(4).partition() != 0)  ITFAILS;

        if ( !soft_equiv(map->get_node(4).point().x(), 1.7) )   ITFAILS;
        if ( !soft_equiv(map->get_node(4).point().y(), 2.7) )   ITFAILS;
        if ( !soft_equiv(map->get_node(4).point().z(), 3.7) )   ITFAILS;
        if ( map->get_node(4).point().handle() != 5)            ITFAILS;
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "One processor map_build test ok on " << node;
        ut.passes( m.str() );
    }
}

//---------------------------------------------------------------------------//
// two processor map building test
void two_pe_map_build_test(Parallel_Unit_Test &ut)
{
    // typedefs
    typedef int                                             OrdinateType;
    typedef double                                          DataType;
    typedef coupler::STAR_File_IO<OrdinateType, DataType>   STAR_File_IO_t;
    typedef typename STAR_File_IO_t::Point_Map_t            Point_Map_t;
    typedef typename Point_Map_t::Point_t                   Point_t; 
    typedef SP<Point_Map_t>                                 SP_Point_Map;
    typedef typename STAR_File_IO_t::Vec_Point              Vec_Point;

    // this test is for two processors
    if (nodes != 2)
        return;

    // make a partitioner
    SP<Partitioner> p;

    // database
    Simple_Partitioner::SP_Std_DB db( new Std_DB("partition") );
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
    }

    // build the mesh
    p->build();

    // make a STAR_File_IO object
    STAR_File_IO_t star_io( get_comm_world(), get_comm_world() );

    // build Star data on node 0
    if (node == 0)
    {
        // make some fake Star data
        Vec_Point points = make_field_pts<OrdinateType, DataType>();

        // assign it to the file_io object
        star_io.set_point_data(points);
    }

    // build the local maps on all nodes
    SP_Point_Map map = star_io.build_map(p);

    // check node 0 local map
    if (node == 0)
    {
        // number of nodes in the map
        if (map->num_nodes() != 3)          ITFAILS;
        
        // unique set of partition ids
        if (map->partitions().size() != 1)  ITFAILS;
        if (map->partitions()[0] != 0)      ITFAILS;

        // point 1
        if (map->get_node(0).partition() != 0)                  ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().x(), 2.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().y(), 2.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().z(), 2.5) )   ITFAILS;
        if (map->get_node(0).point().handle() != 1)             ITFAILS;

        // point 3
        if (map->get_node(1).partition() != 0)                  ITFAILS;
        if ( !soft_equiv(map->get_node(1).point().x(), 2.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(1).point().y(), 7.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(1).point().z(), 2.5) )   ITFAILS;
        if (map->get_node(1).point().handle() != 3)             ITFAILS;

        // point 5
        if (map->get_node(2).partition() != 0 )                 ITFAILS;
        if ( !soft_equiv(map->get_node(2).point().x(), 1.7) )   ITFAILS;
        if ( !soft_equiv(map->get_node(2).point().y(), 2.7) )   ITFAILS;
        if ( !soft_equiv(map->get_node(2).point().z(), 3.7) )   ITFAILS;
        if (map->get_node(2).point().handle() != 5)             ITFAILS;
    }

    // check node 1 local map
    if (node == 1)
    {
        // number of nodes in the map
        if (map->num_nodes() != 2)              ITFAILS;
        
        // unique set of partition ids
        if (map->partitions().size() != 1)      ITFAILS;
        if (map->partitions()[0] != 1)          ITFAILS;

        // point 2
        if (map->get_node(0).partition() != 1)                  ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().x(), 7.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().y(), 2.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().z(), 2.5) )   ITFAILS;
        if (map->get_node(0).point().handle() != 2)             ITFAILS;

        // point 4
        if (map->get_node(1).partition() != 1)                  ITFAILS;
        if ( !soft_equiv(map->get_node(1).point().x(), 7.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(1).point().y(), 7.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(1).point().z(), 2.5) )   ITFAILS;
        if (map->get_node(1).point().handle() != 4)             ITFAILS;
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Two processor map_build test ok on " << node;
        ut.passes(m.str());
    }
}

//---------------------------------------------------------------------------//
// four processor map building test
void four_pe_map_build_test(Parallel_Unit_Test &ut)
{
    // typedefs
    typedef int                                             OrdinateType;
    typedef double                                          DataType;
    typedef coupler::STAR_File_IO<OrdinateType, DataType>   STAR_File_IO_t;
    typedef typename STAR_File_IO_t::Point_Map_t            Point_Map_t;
    typedef typename Point_Map_t::Point_t                   Point_t; 
    typedef SP<Point_Map_t>                                 SP_Point_Map;
    typedef typename STAR_File_IO_t::Vec_Point              Vec_Point;

    // this test is for four processors
    if (nodes != 4)
        return;

    // make a partitioner
    SP<Partitioner> p;

    // database
    Simple_Partitioner::SP_Std_DB db( new Std_DB("partition") );
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
    }

    // build the mesh
    p->build();

    // make a STAR_File_IO object
    STAR_File_IO_t star_io( get_comm_world(), get_comm_world() );

    // build Star data on node 0
    if (node == 0)
    {
        // make some fake Star data
        Vec_Point points = make_field_pts<OrdinateType, DataType>();

        // assign it to the file_io object
        star_io.set_point_data(points);
    }

    // build the local maps on all nodes
    SP_Point_Map map = star_io.build_map(p);

    // check node 0 local map
    if (node == 0)
    {
        // number of nodes in the map
        if (map->num_nodes() != 2)              ITFAILS;
        
        // unique set of partition ids
        if (map->partitions().size() != 1)      ITFAILS;
        if (map->partitions()[0] != 0)          ITFAILS;

        // point 1
        if (map->get_node(0).partition() != 0)                  ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().x(), 2.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().y(), 2.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().z(), 2.5) )   ITFAILS;
        if (map->get_node(0).point().handle() != 1)             ITFAILS;

        // point 5
        if (map->get_node(1).partition() != 0)                  ITFAILS;
        if ( !soft_equiv(map->get_node(1).point().x(), 1.7) )   ITFAILS;
        if ( !soft_equiv(map->get_node(1).point().y(), 2.7) )   ITFAILS;
        if ( !soft_equiv(map->get_node(1).point().z(), 3.7) )   ITFAILS;
        if (map->get_node(1).point().handle() != 5)             ITFAILS;
    }

    // check node 1 local map
    if (node == 1)
    {
        // number of nodes in the map
        if (map->num_nodes() != 1)              ITFAILS;
        
        // unique set of partition ids
        if (map->partitions().size() != 1)      ITFAILS;
        if (map->partitions()[0] != 1)          ITFAILS;

        // point 2
        if (map->get_node(0).partition() != 1)                  ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().x(), 7.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().y(), 2.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().z(), 2.5) )   ITFAILS;
        if (map->get_node(0).point().handle() != 2)             ITFAILS;
    }

    // check node 2 local map
    if (node == 2)
    {
        // number of nodes in the map
        if (map->num_nodes() != 1)              ITFAILS;
        
        // unique set of partition ids
        if (map->partitions().size() != 1)      ITFAILS;
        if (map->partitions()[0] != 2)          ITFAILS;

        // point 3
        if (map->get_node(0).partition() != 2)                  ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().x(), 2.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().y(), 7.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().z(), 2.5) )   ITFAILS;
        if (map->get_node(0).point().handle() != 3)             ITFAILS;
    }

    // check node 3 local map
    if (node == 3)
    {
        // number of nodes in the map
        if (map->num_nodes() != 1)              ITFAILS;
        
        // unique set of partition ids
        if (map->partitions().size() != 1)      ITFAILS;
        if (map->partitions()[0] != 3)          ITFAILS;

        // point 4
        if (map->get_node(0).partition() != 3)                  ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().x(), 7.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().y(), 7.5) )   ITFAILS;
        if ( !soft_equiv(map->get_node(0).point().z(), 2.5) )   ITFAILS;
        if (map->get_node(0).point().handle() != 4)             ITFAILS;
    }

    if (ut.numFails == 0)
    {
        ostringstream m;
        m << "Four processor map build test ok on " << node;
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

        if(nodes == 1)
        {
            int num_blocks_I = 1;
            int num_blocks_J = 1;
            int num_sets     = 1;
            test_geometry_read(num_blocks_I, num_blocks_J, num_sets,
                               "test_star_geom.inp", ut);
            test_power_write(num_blocks_I, num_blocks_J, num_sets,
                             "test_power_write.inp", ut);

            serial_map_build_test(ut);
            ++gpass;
        }
        else if(nodes == 2)
        {
            // 1 Set Tests
            int num_blocks_I = 2;
            int num_blocks_J = 1;
            int num_sets     = 1;
            test_geometry_read(num_blocks_I, num_blocks_J, num_sets,
                               "test_star_geom.inp", ut);
            test_power_write(num_blocks_I, num_blocks_J, num_sets,
                             "test_power_write.inp", ut);
        
            num_blocks_I = 1;
            num_blocks_J = 1;
            num_sets     = 2;
            test_geometry_read(num_blocks_I, num_blocks_J, num_sets,
                               "test_star_geom.inp", ut);
            test_power_write(num_blocks_I, num_blocks_J, num_sets,
                             "test_power_write.inp", ut);

            two_pe_map_build_test(ut);
        }
        else if(nodes == 4)
        {
            // 1 Set Tests
            int num_blocks_I = 2;
            int num_blocks_J = 2;
            int num_sets     = 1;
            test_geometry_read(num_blocks_I, num_blocks_J, num_sets,
                               "test_star_geom.inp", ut);
            test_power_write(num_blocks_I, num_blocks_J, num_sets,
                             "test_power_write.inp", ut);

            // 2 Set Tests
            num_blocks_I = 2;
            num_blocks_J = 1;
            num_sets     = 2;
            test_geometry_read(num_blocks_I, num_blocks_J, num_sets,
                               "test_star_geom.inp", ut);
            test_power_write(num_blocks_I, num_blocks_J, num_sets,
                             "test_power_write.inp", ut);

            // 4 Set Tests
            num_blocks_I = 1;
            num_blocks_J = 1;
            num_sets     = 4;
            test_geometry_read(num_blocks_I, num_blocks_J, num_sets,
                               "test_star_geom.inp", ut);
            test_power_write(num_blocks_I, num_blocks_J, num_sets,
                             "test_power_write.inp", ut);

            four_pe_map_build_test(ut);
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
        std::cout << "ERROR: While testing tstSTAR_Parser, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing tstSTAR_Parser, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstSTAR_Parser.cc
//---------------------------------------------------------------------------//
