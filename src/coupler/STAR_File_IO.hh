//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   coupler/STAR_File_IO.hh
 * \author Gregory Davidson
 * \date   Tue Jun 07 15:15:51 2011
 * \brief  Parses STAR-CCM+ input and output files.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef coupler_STAR_File_IO_hh
#define coupler_STAR_File_IO_hh

#include <string>
#include <vector>

#include "comm/global.hh"
#include "STAR_Data.hh"
#include "Point_Map.hh"

namespace coupler
{

//===========================================================================//
/*!
 * \class STAR_File_IO
 * \brief Parses STAR-CCM+ input and output files.
 */
/*! 
 * \example coupler/test/tstSTAR_File_IO.cc
 *
 * Test of STAR_Parser.
 */
//===========================================================================//
/*
template<class OrdinateType_T, class DataType_T>
class STAR_File_IO
{
  public:
    //@{
    //! Useful typedefs.
    typedef OrdinateType_T                          OrdinateType;
    typedef DataType_T                              DataType;
    typedef denovo::Point<OrdinateType, DataType>   Point_t;
    typedef std::vector<Point_t>                    Vec_Point;
    typedef Point_Map<OrdinateType, DataType>       Point_Map_t;
    typedef denovo::SP<Point_Map_t>                 SP_Point_Map;
    typedef std::vector<Vec_Point>                  Global_Map;
    typedef nemesis::Communicator_t                 Communicator_t;
    typedef denovo::SP<kba::Partitioner>            SP_Partitioner;
    typedef std::string                             FileName;
    typedef std::vector<OrdinateType>               Vec_Ord;
    //@}

  private:
    // Set constant communicator
    const Communicator_t &d_set_const_comm;
    
    // Block constant communicator
    const Communicator_t &d_block_const_comm;

    // Input point data
    Vec_Point d_points;

    // Global Map vector
    Global_Map d_global_map;

  public:
    // Constructor taking just a set constant communicator (assumes 1 set)
    explicit STAR_File_IO(const Communicator_t &set_const_comm);

    // Constructor taking both set and block constant communicators.
    STAR_File_IO(const Communicator_t &set_const_comm, 
                 const Communicator_t &block_const_comm);

    // Read geometry file
    void read_geometry_file(const FileName &filename);

    // Write power file
    void write_power_file(const FileName &filename, 
                          SP_Point_Map    point_map) const;

    // Build local map from geometry file input
    SP_Point_Map build_map(SP_Partitioner partitioner);

    //! Set the point vector externally.
    void set_point_data(const Vec_Point &points)
    {
        d_points = points;
    }

    //! Get the point vector.
    const Vec_Point& get_point_data() const
    {
        return d_points;
    }

  private:
    // >>> PRIVATE TYPEDEFS
    typedef std::vector<char>       Buffer;
    typedef nemesis::Request        Request;

    // >>> IMPLEMENTATION FUNCTIONS
    // Open an input file
    void open_input_file(std::ifstream     &input_stream,
                         const FileName    &filename,
                         const std::string &file_type) const;
    // Open an output file
    void open_output_file(std::ofstream     &output_stream,
                          const FileName    &filename,
                          const std::string &file_type,
                          bool append = false) const;
    // Build the map in serial mode
    SP_Point_Map build_map_serial(SP_Partitioner partitioner);
    // Build the map in parallel mode
    SP_Point_Map build_map_parallel(SP_Partitioner partitioner);
    // Send the sizes of the point vectors to the various nodes
    void send_buffer_sizes(SP_Partitioner partitioner, 
                           Vec_Ord &buffer_sizes);
    // Send the point buffer to the various nodes
    void send_buffers(SP_Partitioner partitioner, 
                      const Vec_Ord &buffer_sizes) const;
    // Fill a map with data from a received buffer
    void fill_map(const Buffer &buffer, SP_Point_Map map) const;
};
*/

} // end namespace coupler


#endif // coupler_STAR_File_IO_hh

//---------------------------------------------------------------------------//
//              end of coupler/STAR_File_IO.hh
//---------------------------------------------------------------------------//
