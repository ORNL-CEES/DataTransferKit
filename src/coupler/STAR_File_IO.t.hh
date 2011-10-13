//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/STAR_File_IO.t.hh
 * \author Gregory Davidson
 * \date   Tue Jun 07 15:15:51 2011
 * \brief  Implements the STAR_File_IO class.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_STAR_File_IO_t_hh
#define coupler_STAR_File_IO_t_hh

#include <fstream>
#include <sstream>
#include <limits>
#include <iomanip>
#include <iterator>

#include "STAR_File_IO.hh"
#include "comm/global.hh"
#include "comm/SpinLock.hh"
#include "utils/Packing_Utils.hh"

namespace coupler
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor taking just a set constant communicator.  This is for
 *        coupling on one set (or with codes that do not use sets).
 */
template<class OrdinateType_T, class DataType_T>
STAR_File_IO<OrdinateType_T, DataType_T>::STAR_File_IO(
    const Communicator_t &set_const_comm)
    : d_set_const_comm(set_const_comm)
    , d_block_const_comm(set_const_comm)
{   }

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor taking both a set constant and block constant
 *        communicator. 
 */
template<class OrdinateType_T, class DataType_T>
STAR_File_IO<OrdinateType_T, DataType_T>::STAR_File_IO(
    const Communicator_t &set_const_comm,
    const Communicator_t &block_const_comm)
    : d_set_const_comm(set_const_comm)
    , d_block_const_comm(block_const_comm)
{   }

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Read in a STAR-CCM+ geometry file.
 */
template<class OrdinateType_T, class DataType_T>
void STAR_File_IO<OrdinateType_T, DataType_T>::read_geometry_file(
    const FileName &filename)
{
    Require ( !filename.empty() );

    // Clear the point vector
    d_points.clear();

    // Only read on Node 0
    if(nemesis::node() == 0)
    {
        // Open file
        std::ifstream input;
        open_input_file(input, filename, "geometry");

        // Skip the first line (header)
        input.ignore(std::numeric_limits<OrdinateType>::max(), '\n');

        // Create some data buffers
        DataType dummy, x, y, z;

        // Create handle
        OrdinateType handle = 0;

        while(!input.eof())
        {
            // Data is stored: Cell_Id, Cell_Index, Region_Index, Volume,
            //                 X, Y, Z.
            // We only need volume, X, Y, and Z.
            input >> dummy >> dummy >> dummy >> dummy >> x >> y >> z;
            
            // Push back the data
            if( input.good() )
            {
                d_points.push_back( Point_t(x, y, z, handle) );
                ++handle;
            }
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write the STAR-CCM+ power file.
 */
template<class OrdinateType_T, class DataType_T>
void STAR_File_IO<OrdinateType_T, DataType_T>::write_power_file(
    const std::string &filename,
    SP_Point_Map       point_map) const
{
    Require (point_map);
    Require ( !filename.empty() );

    // Write the header on Node 0
    if(nemesis::node() == 0)
    {
        std::ofstream output;
        open_output_file(output, filename, "power");
        output << " Enthalpy X Y Z" << std::endl;
    }

    bool do_write = false;
    {
        // >>> We only want the procs on set zero to write.  Figure out if we're
        // on set zero.

        // Set the communicator to a block const communicator.  Now every
        // "node 0" is on set 0.
        nemesis::set_internal_comm(d_block_const_comm);
        if(nemesis::node() == 0)
        {
            // On set zero. Do the write
            do_write = true;
        }
    }

    // Set the communicator to a set const communicator.  Now we can loop over
    // the blocks (using a SpinLock) and write the points.  Only set 0 will
    // write. 
    nemesis::set_internal_comm(d_set_const_comm);

    // Write the data in parallel
    {
        nemesis::HTSyncSpinLock lock;

        if(do_write)
        {
            // Open the file (append mode)
            bool append = true;
            std::ofstream output;
            open_output_file(output, filename, "power", append);

            // Set the precision
            output << std::setprecision(12);

            // Loop over the point map and write to the file
            for(typename Point_Map_t::const_iterator iter = point_map->begin(),
                                                 iter_end = point_map->end();
                iter != iter_end; ++iter)
            {
                output << iter->get("power") << "  " 
                       << iter->point().x() << "  " << iter->point().y() 
                       << "  " << iter->point().z() << std::endl;
            }
        }
    }

    // Reset the communicator
    nemesis::reset_internal_comm();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the local map from input points.
 */
template<class OrdinateType_T, class DataType_T>
typename STAR_File_IO<OrdinateType_T, DataType_T>::SP_Point_Map 
STAR_File_IO<OrdinateType_T, DataType_T>::build_map(SP_Partitioner partitioner)
{
    Require (partitioner);

    if(nemesis::nodes() == 1)
    {
        return build_map_serial(partitioner);
    }
    else
    {
        return build_map_parallel(partitioner);
    }
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Opens file \a filename for input.
 */
template<class OrdinateType_T, class DataType_T>
void STAR_File_IO<OrdinateType_T, DataType_T>::open_input_file(
    std::ifstream     &input_stream,
    const FileName    &filename,
    const std::string &file_type) const
{
    // Open the file
    input_stream.open( filename.c_str() );
    if(!input_stream)
    {
        std::ostringstream error_msg;
        error_msg << "Unable to open STAR-CCM+ " << file_type << " file "
                  << filename << "!";

        Insist( input_stream, error_msg.str() );
    }
}    

//---------------------------------------------------------------------------//
/*!
 * \brief Opens file \a filename for input.
 */
template<class OrdinateType_T, class DataType_T>
void STAR_File_IO<OrdinateType_T, DataType_T>::open_output_file(
    std::ofstream     &output_stream,
    const FileName    &filename,
    const std::string &file_type,
    bool               append) const
{
    // Open the file
    if(append)
    {
        output_stream.open(filename.c_str(), std::ios::out | std::ios::app);
    }
    else
    {
        output_stream.open( filename.c_str() );
    }

    if(!output_stream)
    {
        std::ostringstream error_msg;
        error_msg << "Unable to open STAR-CCM+ " << file_type << " file "
                  << filename << "!";

        Insist( output_stream, error_msg.str() );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the local map from input points in serial mode.
 */
template<class OrdinateType_T, class DataType_T>
typename STAR_File_IO<OrdinateType_T, DataType_T>::SP_Point_Map
STAR_File_IO<OrdinateType_T, DataType_T>::build_map_serial(
    SP_Partitioner partitioner)
{
    Require (partitioner);

    // make a point map on this node
    SP_Point_Map map( new Point_Map_t(partitioner) );
    Ensure (map);

    // add each point to the map
    for (typename Vec_Point::const_iterator point_iter = d_points.begin(),
                                        point_iter_end = d_points.end();
         point_iter != point_iter_end; ++point_iter)
    {
        map->add_point(*point_iter);
    }

    // complete the mapping
    map->complete();

    // return the map
    return map;
}


//---------------------------------------------------------------------------//
/*!
 * \brief Build the local map from input points in parallel mode.
 */
template<class OrdinateType_T, class DataType_T>
typename STAR_File_IO<OrdinateType_T, DataType_T>::SP_Point_Map 
STAR_File_IO<OrdinateType_T, DataType_T>::build_map_parallel(
    SP_Partitioner partitioner)
{
    typedef nemesis::Request            Request;
    typedef std::vector<char>           Buffer;

    Require(partitioner);
    Check ( nemesis::nodes() == 
            partitioner->num_sets() * partitioner->num_blocks() );

    // make a point map on this node
    SP_Point_Map map( new Point_Map_t(partitioner) );
    
    // Get the node id
    OrdinateType node = nemesis::node();

    // >>> REQUEST THE BUFFER SIZE
    OrdinateType this_buffer_size;
    Request size_request = nemesis::receive_async(&this_buffer_size, 1, 0);

    // >>> SEND THE BUFFER SIZES OUT FROM NODE 0
    Vec_Ord buffer_sizes( partitioner->num_blocks() );
    if(node == 0)
    {
        send_buffer_sizes(partitioner, buffer_sizes);
    }

    // >>> WAIT FOR THE BUFFER SIZES TO BE RECEIVED
    size_request.wait();

    // >>> POST RECEIVE OF BUFFER    
    // Make the buffer
    Buffer buffer(this_buffer_size);
    Request buffer_request;
    if(this_buffer_size > 0)
    {
        // Receive the buffer
        buffer_request = nemesis::receive_async(&buffer[0], this_buffer_size, 0);
    }

    // >>> SEND THE BUFFERS
    if(node == 0)
    {
        send_buffers(partitioner, buffer_sizes);
    }

    if(this_buffer_size > 0)
    {
        // >>> WAIT UNTIL THE BUFFER COMES IN
        buffer_request.wait();

        // >>> FILL THE MAP
        fill_map(buffer, map);
    }

    // complete the mapping
    map->complete();

    // Communication is done.  Do barrier
    nemesis::global_barrier();

    // return the map
    return map;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send the map buffer sizes to the other nodes.
 *
 * \note This function should be executed only by node 0.
 */
template<class OrdinateType_T, class DataType_T>
void STAR_File_IO<OrdinateType_T, DataType_T>::send_buffer_sizes(
    SP_Partitioner  partitioner,
    Vec_Ord        &buffer_sizes)
{
    typedef std::vector<char>               Buffer;
    
    Require (partitioner);

    // set the size of the global map vector
    d_global_map.resize( partitioner->num_blocks() );

    // add each point to the global map
    for (typename Vec_Point::const_iterator point_iter = d_points.begin(),
                                        point_iter_end = d_points.end();
         point_iter != point_iter_end; ++point_iter)
    {
        OrdinateType block = 
            partitioner->point_query( point_iter->coords() );

        // we require that the point is inside the global mesh
        Check(block != -1);
         
        d_global_map[block].push_back(*point_iter);
    }
        
    // compute the size of the star data buffers going to each block
    for (typename Global_Map::const_iterator iter = d_global_map.begin(),
                                         iter_end = d_global_map.end(); 
         iter != iter_end; ++iter)
    {
        // Calculate what block we are sending to
        OrdinateType block = iter - d_global_map.begin();
        Check ( block < buffer_sizes.size() );

        // Create a packer
        denovo::Packer p;

        // Pack in point vector to compute the size
        p.compute_buffer_size_mode();
        for (typename Vec_Point::const_iterator vec_iter = iter->begin(),
                                            vec_iter_end = iter->end();
             vec_iter != vec_iter_end; ++vec_iter)
        {
            p << *vec_iter;
        }
        buffer_sizes[block] = p.size();

        // Loop over the sets and send the points to all the nodes in the set
        for(OrdinateType set = 0, set_end = partitioner->num_sets(); 
            set < set_end; ++set)
        {
            // Get the destination node
            OrdinateType destination = partitioner->node_query(block, set);
            Check(destination >= 0);
            Check( destination < nemesis::nodes() );

            // Send the buffer size
            nemesis::send_async(&buffer_sizes[block], 1, destination);
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send the map information to the other nodes.
 *
 * \note This function should be executed only by node 0.
 */
template<class OrdinateType_T, class DataType_T>
void STAR_File_IO<OrdinateType_T, DataType_T>::send_buffers(
    SP_Partitioner partitioner, const Vec_Ord& buffer_sizes) const
{
    Require (partitioner);

    // compute the size of the star data buffers going to each block
    for (typename Global_Map::const_iterator iter = d_global_map.begin(),
                                         iter_end = d_global_map.end(); 
         iter != iter_end; ++iter)
    {
        // Calculate what block we are sending to
        OrdinateType block = iter - d_global_map.begin();
        Check ( block < buffer_sizes.size() );

        // create a packer
        denovo::Packer p;

        // pack the point vector into the buffer
        Buffer buffer(buffer_sizes[block]);
        if(buffer_sizes[block] > 0)
        {
            p.set_buffer(buffer_sizes[block], &buffer[0]);
            for (typename Vec_Point::const_iterator vec_iter = iter->begin(),
                                                vec_iter_end = iter->end();
                 vec_iter != vec_iter_end; ++vec_iter)
            {
                p << *vec_iter;
            }
        }

        // Loop over the sets and send the points to all the nodes in the set
        for(OrdinateType set = 0, set_end = partitioner->num_sets(); 
            set < set_end; ++set)
        {
            // Get the destination node
            OrdinateType destination = partitioner->node_query(block, set);
            Check(destination >= 0);
            Check( destination < nemesis::nodes() );

            // Send the buffer
            nemesis::send_async(&buffer[0], buffer_sizes[block], destination);
        }
    }
}
    
//---------------------------------------------------------------------------//
/*!
 * \brief Fill the map with information from a received buffer.
 */
template<class OrdinateType_T, class DataType_T>
void STAR_File_IO<OrdinateType_T, DataType_T>::fill_map(
    const Buffer &buffer, SP_Point_Map map) const
{
    Require (map);

    // compute the number of points in the buffer
    OrdinateType num_points = buffer.size() / sizeof(Point_t);

    // unpack the buffer into star data containers and build the 
    // local map
    denovo::Unpacker u;
    u.set_buffer(buffer.size(), &buffer[0]);
 
    for (OrdinateType l = 0; l < num_points; ++l)
    {
        // Make a dummy point
        Point_t dummy_point(0.0, 0.0, 0.0, 0);

        // Unpack into this point
        u >> dummy_point;

        // Add the point to the map
        map->add_point(dummy_point);
    }
}

} // end namespace coupler


#endif // coupler_STAR_File_IO_t_hh


//---------------------------------------------------------------------------//
//                 end of STAR_File_IO.t.hh
//---------------------------------------------------------------------------//
