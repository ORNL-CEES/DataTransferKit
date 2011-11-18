//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Mapper.hh
 * \author Stuart Slattery
 * \date   Tue Nov 08 12:31:05 2011
 * \brief  Mapper class defintion.
 * \note   Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.hh,v 1.4 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef core_Mapper_hh
#define core_Mapper_hh

#include <string>
#include <list>
#include <vector>

#include "Transfer_Data_Field.hh"
#include "Transfer_Map.hh"
#include "Message_Buffer.hh"
#include "LG_Indexer.hh"
#include "utils/Packing_Utils.hh"
#include "comm/global.hh"

#include "Teuchos_RCP.hpp"

namespace coupler
{

//===========================================================================//
/*!
 * \class Mapper
 * \brief A mapper class to generate a one way mapping for transfer between a
 * source physics and a target physics.
 *
 * \sa Mapper.cc for detailed descriptions.
 */
/*! 
 * \example core/test/tstMapper.cc
 *
 * Test of Mapper.
 */
//===========================================================================//

template<class DataType_T>
class Mapper 
{
  public:

    //@{
    //! Useful typedefs.
    typedef DataType_T                            DataType;
    typedef int                                   HandleType;
    typedef int                                   OrdinateType;
    typedef double                                CoordinateType;
    typedef Message_Buffer<OrdinateType>          Message_Buffer_t;
    typedef typename Message_Buffer_t::Buffer     Buffer;
    typedef std::list<Message_Buffer_t>           BufferList;
    typedef typename BufferList::iterator         BufferList_Iterator;
    typedef Transfer_Data_Field<DataType>         Transfer_Data_Field_t;
    typedef Teuchos::RCP<Transfer_Data_Field_t>   RCP_Transfer_Data_Field;
    typedef nemesis::Communicator_t               Communicator;
    typedef Teuchos::RCP<Transfer_Map>            RCP_Transfer_Map;
    typedef typename Transfer_Map::Map_Iterator   Map_Iterator;
    typedef typename Transfer_Map::Map_Pair       Map_Pair;
    typedef typename Transfer_Map::Set_Iterator   Set_Iterator;
    typedef typename Transfer_Map::Set_Pair       Set_Pair;
    //@}

  public:

    // Constructor.
    Mapper();

    // Destructor.
    ~Mapper();

    // Map the field from the source onto the target.
    void map(const Communicator &comm_global,
	     RCP_Transfer_Data_Field transfer_data_field,
	     RCP_Transfer_Map transfer_map);

  private:

    // Source physics post receives for buffer sizes.
    void source_post_receive_size(LG_Indexer &target_indexer,
				  BufferList &buffer_size_list);

    // Target physics sends point sizes to source.
    void target_send_point_size(RCP_Transfer_Data_Field transfer_data_field,
				LG_Indexer &source_indexer,
				std::vector<CoordinateType> &coordinates,
				std::vector<HandleType> &handles);

    // Source physics process requests for message sizes and post receives
    // for buffers.
    void source_post_receive_buffer(LG_Indexer &target_indexer,
				    BufferList &buffer_size_list,
				    BufferList &buffer_list);

    // Target send points to source.
    void target_send_points(LG_Indexer &source_indexer,
			    const std::vector<CoordinateType> &coordinates,
			    const std::vector<HandleType> &handles);

    // Source physics process request and build part of the map.
    void source_process_points(RCP_Transfer_Data_Field transfer_data_field, 
			       BufferList &buffer_list,
			       RCP_Transfer_Map transfer_map);

    // Target physics post receives for return buffer size.
    void target_post_receive_size(LG_Indexer &source_indexer,
				  BufferList &buffer_size_list);

    // Source physics sends back the number of points it found in its domain
    // back to the target.
    void source_send_point_size(LG_Indexer &target_indexer,
				RCP_Transfer_Map transfer_map);

    // Target physics process request for message sizes and post receives.
    void target_post_receive_buffer(BufferList &buffer_size_list,
				    BufferList &buffer_list);

    // Source physics sends its point handles to the targets.
    void source_send_handles(RCP_Transfer_Map transfer_map);
    
    // Target physics processes handle requests and completes the mapping.
    void target_process_handles(BufferList &buffer_list,
				RCP_Transfer_Map transfer_map);
};

} // end namespace coupler

#endif // core_Mapper_hh

//---------------------------------------------------------------------------//
//              end of core/Mapper.hh
//---------------------------------------------------------------------------//
