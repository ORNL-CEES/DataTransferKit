//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/Coupler.t.hh
 * \author Stuart R. Slattery
 * \date   Fri Jun 10 08:56:27 2011
 * \brief  Coupler class template member definitions
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//


#ifndef coupler_Coupler_t_hh
#define coupler_Coupler_t_hh

#include "Coupler.hh"
#include "STAR_File_IO.hh"
#include "STAR_Data.hh"
#include "Map_Node.hh"
#include "Point_Map.hh"
#include "Messenger.hh"
#include "fields/Field_View.hh"
#include "comm/global.hh"
#include "harness/DBC.hh"
#include "utils/SP.hh"
#include "utils/Packing_Utils.hh"
#include "utils/Vector_Lite.hh"

#include <string>
#include <vector>
#include <algorithm>

namespace coupler
{

//===========================================================================//
// GENERAL COUPLER IMPLEMENTATION
//===========================================================================//


//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Coupler constructor.
 *
 * \param comm_world         Communicator that encompasses all of the master
 *                           and the slave application being coupled.
 * \param comm_master        Communicator for the master application.
 * \param comm_slave         Communicator for the slave application.
 * \param external_app_ptr   Smart pointer to the master application.
 * \param neutronics_app_ptr Smart pointer to the slave application.
 */
template<class Master_App, class Slave_App>
Coupler<Master_App, Slave_App>::Coupler(
    Communicator_t comm_world, Communicator_t comm_master, 
    Communicator_t comm_slave, denovo::SP<Master_App> master_app_ptr,
    denovo::SP<Slave_App> slave_app_ptr)
    : d_world_comm(comm_world)
    , d_master_comm(comm_master)
    , d_slave_comm(comm_slave)
    , d_master_app(master_app_ptr)
    , d_slave_app(slave_app_ptr)
{
    nemesis::set_internal_comm(d_comm_world);

    if(d_master_app || d_slave_app)
    {
        d_indexer_master = new LG_Indexer_t(d_world_comm, d_master_comm,
                                            d_master_app);
        Ensure (d_indexer_master);

        d_indexer_slave = new LG_Indexer_t(d_world_comm, d_slave_comm,
                                           d_slave_app);
        Ensure (d_indexer_slave);
    }
    // reset the internal communicator
    nemesis::reset_internal_comm();
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Register a vector of interleaved point coordinates with the Coupler
 *        and their associated handles and sort the resulting Point objects by 
 *        handle. 
 *
 * \param points  A vector of doubles containing interleaved point coordinates.
 *                (i.e. \f$ [x_1, y_1, z_1, x_2, y_2, z_2, ... , x_N, y_N, 
 *                           z_N] \f$ ) 
 * \param handles A vector of integer handles corresponding to the points.
 */
void Coupler::register_points(const Vec_DataType &points, 
                              const Vec_Ordinates &handles)
{
    Require( points.size() % 3 == 0 );
    Require( points.size() / 3 == handles.size() );

    OrdinateType num_points = points.size() / 3;
    d_points.reserve(num_points);

    // Add the points the point vector
    for (OrdinateType i = 0; i < num_points; ++i)
    {
        Point_Dbl pt(points[3*i], points[3*i+1], points[3*i+2], handles[i]);
        d_points.push_back(pt);
    }

    // Sort the points by handle
    std::sort(d_points.begin(), d_points.end(), &handle_compare);

    Ensure( d_points.size() == handles.size() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a local map on all processors.
 *
 * This function is for general parallel-parallel transfer. It builds the
 * local map and generates the messenger object on both neutronics and the
 * external application.
 */
void Coupler::build_map()
{
    // operate on the ordered union communicator
    nemesis::set_internal_comm(d_communicator);

    // initialize the local maps
    if (d_master_app)
    {
        d_master_map = new Point_Map();
        Ensure (d_external_map);
    }

    if (d_neutronics_app)
    {
        d_neutronics_map = new Point_Map();
        Ensure (d_neutronics_map);
    }

    // set local process id iteration bounds
    int neutronics_begin, neutronics_end;
    int external_begin, external_end;
    if (d_neutronics_app || d_external_app)
    {
        neutronics_begin = 0;
        neutronics_end = d_indexer_neutronics->size();
        external_begin = 0;
        external_end = d_indexer_external->size();
    }

    // send points from the external application to neutronics for neutronics
    // local map generation
    if (d_external_app)
    {
        // create a packer to send points to neutronics
        denovo::Packer p;

        // compute the size of the buffer
        p.compute_buffer_size_mode();
        for (Vec_Point::const_iterator iter = d_points.begin(),
                                   iter_end = d_points.end();
             iter != iter_end; ++iter)
        {
            p << *iter;
        }
        int buffer_size = p.size();

        // pack the points into the buffer
        Buffer buffer(buffer_size);
        if (buffer_size > 0)
        {
            p.set_buffer(buffer_size, &buffer[0]);
            for (Vec_Point::const_iterator iter = d_points.begin(),
                                       iter_end = d_points.end();
                 iter != iter_end; ++iter)
            {
                p << *iter;
            }
        }
        Check( buffer.size() == buffer_size );

        // non-blocking send the local points to all neutronics processes
        for (int i = neutronics_begin; i < neutronics_end; ++i)
        {
            // get the global index for neutronics that we are sending to
            int destination = d_indexer_neutronics->l2g(i);
	    
            // send a message so neutronics knows what is coming
            nemesis::send_async(&buffer_size, 1, destination);

            // send the actual buffer
            if (buffer_size > 0)
            {
                nemesis::send_async(&buffer[0], buffer_size, destination);
            }
        }
    }

    // neutronics receives the points and builds a local map
    if (d_neutronics_app)
    {
        // every external process sends neutronics it's points
        for (int i = external_begin; i < external_end; ++i)
        {
            // get the global index for external that we are receiving from
            int source = d_indexer_external->l2g(i);

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
                    if ( d_neutronics->partitioner()->get_mesh()->find_cell(
                             point.get_coords(), dummy_vector) )
                    {
                        // Make a new map node
                        Map_Node map_node(point.get_handle(), source, point);

                        // add it to the neutronics map
                        d_neutronics_map->add_map_node(map_node);
                    }
                }
            }
        }
	    
        // complete the neutronics mapping
        d_neutronics_map->complete();
    }

    // barrier after the neutronics map generation
    nemesis::global_barrier();

    // send points back to the external application from neutronics to
    // complete the local map on the external application
    if (d_neutronics_app)
    {
        // get the unique partition id vector for this neutronics node
        Vec_Int partitions = d_neutronics_map->get_partitions();

        // get the number of points this neutronics node is sending to each
        // external process
        Vec_Int num_send_points = d_neutronics_map->get_partition_pts();

        Check( partitions.size() == num_send_points.size() );

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
                p << d_neutronics_map->get_all_nodes()[node_counter].
                    get_node_point();
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
                p << d_neutronics_map->get_all_nodes()[node_counter].
                    get_node_point();
                ++node_counter;
            }

            // get the global index to external that we are sending to
            int destination = partitions[i];

            // send the size of the buffer
            nemesis::send_async(&buffer_size, 1, destination);

            // send the buffer
            nemesis::send_async(&buffer[0], buffer_size, destination);
        }

        // make sure we sent all of the node points
        Check( node_counter == d_neutronics_map->get_num_nodes() );
    }

    // the external application gets points from neutronics to complete its map
    if (d_external_app)
    {
        // get a buffer from all neutronics processes
        for (int i = neutronics_begin; i < neutronics_end; ++i)
        {
            // get the global index to neutronics that we are receiving from
            int source = d_indexer_neutronics->l2g(i);
	    
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
                    Vec_Point::iterator it = 
                        std::lower_bound(d_points.begin(), d_points.end(),
                                         dummy, &handle_compare);

                    // add a node to the map for the point
                    Map_Node map_node(it->get_handle(), source, *it);
                    d_external_map->add_map_node(map_node);
                }
            }
        }

        // complete the external application mapping
        d_external_map->complete();
    }

    // barrier after the external map generation
    nemesis::global_barrier();

    // create a messenger on each process
    if (d_neutronics_app)
    {
        d_messenger_neut = new Messenger(d_communicator, d_neutronics_map);
        Ensure (d_messenger_neut);

        d_messenger_neut->setup();
    }
    if (d_external_app)
    {
        d_messenger_ext = new Messenger(d_communicator, d_external_map);
        Ensure (d_messenger_ext);

        d_messenger_ext->setup();
    }

    // reset the internal communicator
    nemesis::reset_internal_comm();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build a local map on all processors.
 *
 * This function is overloaded for a STAR-CCM+ geometry input file.
 *
 * \param STAR_geom_file STAR-CCM+ formatted text file containing geometric
 *                       information.
 */
void Coupler::build_map(const std::string &STAR_geom_file)
{
    // make a STAR_File_IO object
    STAR_File_IO star_io(d_neutronics_set_comm, d_neutronics_block_comm);

    // read the STAR geometry file
    star_io.read_geometry_file(STAR_geom_file);

    // build the local map
    d_neutronics_map = star_io.build_map(d_neutronics->partitioner());

    // populate the local vector of points
    const std::vector<Map_Node> &map_nodes = d_neutronics_map->get_all_nodes();
    for (std::vector<Map_Node>::const_iterator iter = map_nodes.begin(),
                                           iter_end = map_nodes.end(); 
         iter != iter_end; ++iter)
    {
        d_points.push_back( iter->get_node_point() );
        d_volumes.push_back( iter->get_node_volume() );
    }

    Ensure( d_points.size() == d_neutronics_map->get_num_nodes() );
    Ensure( d_volumes.size() == d_neutronics_map->get_num_nodes() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Transfer power data from Neutronics to an external application.
 *
 * This function is overloaded for general many to many parallel transfer.
 */
void Coupler::transfer_power()
{
    // neutronics sends the power
    if (d_neutronics_app)
    {
        // make a vector of powers to populate
        Vec_Dbl powers(d_neutronics_map->get_num_nodes(), 0.0);
        View_Field power_view(powers);

        // make a dummy vector of volumes - should remove after API change for
        // the power calculator
        d_volumes.resize(d_neutronics_map->get_num_nodes());
        std::fill(d_volumes.begin(), d_volumes.end(), 1.0);

        // make a vector of the local neutronics points in the map
        Vec_Point neutronics_points;
        for (int i = 0; i < powers.size(); ++i)
        {
            neutronics_points.push_back(
                d_neutronics_map->get_all_nodes()[i].get_node_point() );
        }
        Check( neutronics_points.size() == d_neutronics_map->get_num_nodes() );
	      
        // compute the power at all local points     
        d_neutronics->calc_power(neutronics_points, d_volumes, power_view);


        // >>>>>>>>>> Temporary
        for(View_Field::const_iterator iter = power_view.begin();
            iter != power_view.end(); ++iter)
        {
            std::cout << *iter << std::endl;
        }
        // >>>>>>>>>> Temporary


        // assign the power values to the map nodes
        for (int j = 0; j < powers.size(); ++j)
        {
            d_neutronics_map->get_all_nodes()[j].set_node_power( powers[j] );
        }

        // send the power
        d_messenger_neut->send_power();
    }

    // external application receives the power
    if (d_external_app)
    {
        d_messenger_ext->receive_power();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Transfer power data from Neutronics to an external application.
 *
 * This function is overloaded for coupling to STAR-CCM+.
 *
 * \param STAR_power_file STAR-CCM+ formatted output file with power data.
 */
void Coupler::transfer_power(const std::string &STAR_power_file)
{
    // create a power vector field view
    std::vector<double> powers(d_points.size(), 0.0);
    View_Field power_view(powers);

    // compute the power at all points
    d_neutronics->calc_power(d_points, d_volumes, power_view);

    // create a STAR_File_IO object
    STAR_File_IO star_io(d_neutronics_set_comm, d_neutronics_block_comm);

    // write out a file with the powers
    star_io.write_power_file(STAR_power_file, d_points, powers);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the vector of powers registered with the Coupler.
 *
 * The returned vector is sorted in increasing order of the point handle that
 * is associated with the power value.
 *
 * \param powers Reference to a vector<double> that will be populated with
 *               power data.
 */
void Coupler::get_power(Vec_Dbl &powers)
{
    if (d_external_app)
    {
        Require( powers.size() == d_external_map->get_num_nodes() );

        // sort the map nodes in the external map by handle
        d_external_map->sort_nodes_by_handle();

        // fill the power vector with the powers in the map nodes
        for (int i = 0; i < d_external_map->get_num_nodes(); ++i)
        {
            powers[i] = d_external_map->get_all_nodes()[i].get_node_power();
        }

        // resort the map nodes in the external map by partition for future data
        // transfers 
        d_external_map->sort_nodes_by_partition();
    }
} 


} // end namespace coupler


#endif // coupler_Coupler_t_hh

//---------------------------------------------------------------------------//
//                 end of Coupler.t.hh
//---------------------------------------------------------------------------//
