//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Mapper.cc
 * \author Stuart Slattery
 * \date   Tue Nov 08 12:31:05 2011
 * \brief  Mapper class member definitions.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include <algorithm>

#include "Mapper.hh"
#include "harness/DBC.hh"

namespace coupler
{
//---------------------------------------------------------------------------//
// CONSTRUCTOR and DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Mapper::Mapper(const Communicator &comm_global,
	       const std::string &field_name,
	       SP_Physics source,
	       SP_Physics target)
    : d_comm_global(comm_global)
    , d_field_name(field_name)
    , d_source(source)
    , d_target(target)
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Mapper::~Mapper()
{ /* ... */ }


//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Map the field from the source onto the target.
 */
void Mapper::map()
{
    //  Set the internal communicator.
    nemesis::set_internal_communicator(d_comm_global);

    // Initialize a map.
    SP_Transfer_Map new_map = new Transfer_Map();

    // Create an empty list of message buffers.
    BufferList buffer_size_list;
    BufferList buffer_list;

    // Source physics post receives for the message buffer size.
    if ( d_source->te() )
    {
	source_post_receive_size(new_map, buffer_size_list);
    }

    // Target point coordinate vector iterators.
    Coord_Iterator points_begin, points_end;

    // Target point handle vector iterators.
    Handle_Iterator handles_begin, handles_end;

    // Target physics sends message buffer sizes to source.
    if ( d_target->te() )
    {
	target_send_point_size(points_begin, points_end,
			       handles_begin, handles_end);
    }

    // Source physics processes requests and posts receives for the buffers.
    if ( d_source->te() )
    {
	source_post_receive_buffer(buffer_size_list, buffer_list);
    }

    // Target send buffers with points to the source.
    if ( d_target->te() )
    {
	target_send_points(points_begin, points_end,
			   handles_begin, handles_end);
    }

    // Source physics processes requests of target points from the target
    // physics and builds the topology map. 
    if (d_source->te())
    {
	source_process_points(buffer_list, new_map);
    }

    // Barrier after sending target points from the target physics to the
    // source physics.
    nemesis::global_barrier();

    // Clear the buffer lists to be safe.
    buffer_size_list.clear();
    buffer_list.clear();

    // Target physics post receives for the return communication of which
    // handles were found in the source domain. These are requests are for the
    // size of that communication.
    if ( d_target->te() )
    {
	target_post_receive_size(buffer_size_list);
    }

    // Source physics sends the number of points it found in its domain back
    // to the target.
    if ( d_source->te() )
    {
	source_send_point_size();
    }

    // Target physics processes requests for the number of points from the
    // source and then posts receives if that number is greater than zero.
    if ( d_target->te() )
    {
	target_post_receive_buffer(buffer_size_list, buffer_list);
    }

    // Source physics sends target point handles it found in its domain back
    // to the target.
    if ( d_source->te() )
    {
	source_send_handles(new_map);
    }

    // Target physics processes requests for handles from the source and
    // completes the map.
    if (d_target->te())
    {
	target_process_handles(buffer_list, new_map);
    }

    // Barrier after completion.
    nemesis::global_barrier();

    // Add the new map to the source physics database.
    d_source->set_map(target_physics, field_name, new_map);
    
    // Reset the internal communicator.
    nemesis::reset_internal_comm();
}


//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
// Source physics post receives for buffer sizes.
void Mapper::source_post_receive_size(BufferList &buffer_size_list)
{
    // Initialize.
    Message_Buffer_t &buffer;
    int buffer_size = 0;

    // Post receives for each target process.
    OrdinateType src;
    OrdinateType begin_target = 0;
    OrdinateType end_target = d_target->indexer()->size();
    
    for ( src = begin_target; src != end_target; ++src)
    {
	Check ( src < nemesis::nodes() );

	// Create the buffer and add it to the list.
	buffer_size_list.push_back( Message_Buffer_t(src, 1) );

	// Get the request buffer.
	buffer = buffer_size_list.back();

	// Post asynchronous receive with this receive buffer.
	nemesis::receive_async(buffer.request(),
			       &buffer.buffer()[0],
			       buffer.buffer().size(),
			       buffer.ordinate());
    }

    // Make sure we made all of the buffers we're going to receive.
    Check ( buffer_size_list.size() == d_target->indexer()->size() );
}

//---------------------------------------------------------------------------//
// Target physics sends point sizes to source.
void Mapper::target_send_point_size(Coord_Iterator &points_begin,
			    Coord_Iterator &points_end,
			    Handle_Iterator &handles_begin,
			    Handle_Iterator &handles_end)
{
    // Target physics registers its target points for the field being mapped.
    d_target->te()->register_xyz(field_name, 
				 points_begin, points_end, 
				 handles_begin, handles_end);

    Check( std::distance(points_begin, points_end) % 3 == 0 );
    Check( std::distance(points_begin, points_end) / 3 == 
	   std::distance(handles_begin, handles_end) );

    // Create a packer.
    denovo::Packer p;

    // Compute the size of the buffer.
    Coord_Iterator coord_iter = points_begin;
    Handle_Iterator handle_iter;
    p.compute_buffer_size_mode();

    for (handle_iter = handles_begin, handle_iter != handles_end; ++iter)
    {
	p << *coord_iter;
	coord_iter++;
	p << *coord_iter;
	coord_iter++;
	p << *coord_iter;
	coord_iter++;
	p << *handle_iter;
    }

    int buffer_size = p.size();

    // Send the local target points to all processes of the source physics.
    OrdinateType destination;
    OrdinateType begin_source = 0;
    OrdinateType end_source = d_source->indexer()->size();

    for (int i = begin_source; i < end_source; ++i)
    {
	// Get the global index for the source physics that the buffer is
	// going to.
	destination = d_source->indexer()->l2g(i);

	// Send a message to the source with the size of the buffer that
	// it will get. 
	nemesis::send_async(&buffer_size, 1, destination);
    }
}

//---------------------------------------------------------------------------//
// Source physics process requests for message sizes and post receives
// for buffers.
void Mapper::source_post_receive_buffer(BufferList &buffer_size_list,
				BufferList &buffer_list)
{
    // Initialize.
    OrdinateType src;
    BufferList_Iterator buffer_iter;
    Buffer &buffer;
    int buffer_size;
    denovo::Unpacker u;
	
    // Keep going until all requests have been processed.
    while( !buffer_list.empty() )
    {
	// Find a buffer with a completed communication request.
	buffer_iter = std::find_if(buffer_size_list.begin(), 
				   buffer_size_list.end(), 
				   &Message_Buffer_t::complete);

	// If a completed communication request was found, process it.
	if( buffer_iter != buffer_size_list.end() )
	{
	    // Get the source partition.
	    src = buffer_iter->ordinate();

	    // Get the buffer.
	    buffer = buffer_iter->buffer();

	    // Set the buffer for the unpacker.
	    u.set_buffer( buffer.size(), &buffer[0] );

	    // Get the size of the next buffer we will receive.
	    u >> buffer_size;

	    // Remove this buffer from the size list.
	    buffer_size_list.erase(buffer_iter);

	    // Create the buffer and add it to the list.
	    buffer_list.push_back( Message_Buffer_t(src, buffer_size) );

	    // Clear the buffer just to be safe.
	    buffer.clear();

	    // Get the request buffer.
	    buffer = buffer_list.back();

	    // Post asynchronous receive with this receive buffer.
	    nemesis::receive_async(buffer.request(),
				   &buffer.buffer()[0],
				   buffer.buffer().size(),
				   buffer.ordinate());
	}

	// Make sure we made all of the buffers we're going to receive.
	Check ( buffer_list.size() == d_target->indexer()->size() );
    }
}

//---------------------------------------------------------------------------//
// Target send points to source.
void Mapper::target_send_points(Coord_Iterator points_begin,
				Coord_Iterator points_end,
				Handle_Iterator handles_begin,
				Handle_Iterator handles_end)
{
    // Build a buffer of the local points to send to the source physics.
    Buffer buffer;
    
    // Compute the size of the buffer.
    Coord_Iterator coord_iter = points_begin;
    Handle_Iterator handle_iter;
    p.compute_buffer_size_mode();
    for (handle_iter = handles_begin, handle_iter != handles_end; ++iter)
    {
	p << *coord_iter;
	coord_iter++;
	p << *coord_iter;
	coord_iter++;
	p << *coord_iter;
	coord_iter++;
	p << *handle_iter;
    }
    int buffer_size = p.size();

    // Set the size of the buffer.
    buffer.resize(buffer_size);

    // Pack the data into the buffer.
    Check( buffer.size() == buffer_size );
    if (buffer_size > 0)
    {
	p.set_buffer(buffer_size, &buffer[0]);
	for (handle_iter = handles_begin;
	     handle_iter != handles_end; 
	     ++handle_iter)
	{
	    p << *coord_iter;
	    coord_iter++;
	    p << *coord_iter;
	    coord_iter++;
	    p << *coord_iter;
	    coord_iter++;
	    p << *handle_iter;
	}
    }
    int buffer_size = buffer.size();

    // Send the local target points to all processes of the source physics.
    int destination;
    for (int i = begin_source; i < end_source; ++i)
    {
	// Get the global index for the source physics that the buffer is
	// going to.
	destination = d_source->indexer()->l2g(i);

	// Send the buffer of points.
	if (buffer_size > 0)
	{
	    nemesis::send_async(&buffer[0], buffer_size, destination);
	}
    }
}

//---------------------------------------------------------------------------//
// Source physics process request and build part of the map.
void Mapper::source_process_points(BufferList &buffer_list,
				   SP_Transfer_Map new_map)
{
    // Initialize.
    OrdinateType src;
    BufferList_Iterator buffer_iter;
    Buffer &buffer;
    int buffer_size;
    int num_points;
    int j;
    denovo::Unpacker u;
    HandleType handle;
    Coordinate x, y, z;

    while ( buffer_list.empty() )
    {
	// Find a buffer with a completed communication request.
	buffer_iter = std::find_if(buffer_list.begin(), 
				   buffer_list.end(), 
				   &Message_Buffer_t::complete);

	// If a completed communication request was found, process it.
	if( buffer_iter != buffer_size_list.end() )
	{
	    // Get the source partition.
	    src = buffer_iter->ordinate();

	    // Get the buffer.
	    buffer = buffer_iter->buffer();

	    // Get the buffer size.
	    buffer_size = buffer.size();

	    // Unpack the points and add them to the map.
	    if (buffer_size > 0)
	    {
		// Compute the number of points in the buffer.
		num_points = buffer_size / ( sizeof(double) * 3 + sizeof(int));

		// Unpack the buffer.
		u.set_buffer(buffer_size, &buffer[0]);

		for (j = 0; j < num_points; ++j)
		{
		    u >> x;
		    u >> y;
		    u >> z;
		    u >> handle;

		    // See if this point is in the spatial domain of this
		    // source process.
		    if ( d_source->te()->find_xyz(x, y, z, handle) )
		    {
			// Add the handle to the map with the target rank.
			new_map->add_domain_pair(src, handle);
		    }
		}

		// Remove this buffer from the list.
		buffer_list.erase(buffer_iter);
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// Target physics post receives for return buffer size.
void Mapper::target_post_receive_size(BufferList &buffer_size_list)
{
    // Initialize.
    Message_Buffer_t &buffer;
    int buffer_size = 0;

    // Post receives for each source process.
    OrdinateType src;
    OrdinateType begin_source = 0;
    OrdinateType end_source = d_source->indexer()->size();
    
    for ( src = begin_source; src != end_source; ++src)
    {
	Check ( src < nemesis::nodes() );

	// Create the buffer and add it to the list.
	buffer_size_list.push_back( Message_Buffer_t(src, 1) );

	// Get the request buffer.
	buffer = buffer_size_list.back();

	// Post asynchronous receive with this receive buffer.
	nemesis::receive_async(buffer.request(),
			       &buffer.buffer()[0],
			       buffer.buffer().size(),
			       buffer.ordinate());
    }

    // Make sure we made all of the buffers we're going to receive.
    Check ( buffer_size_list.size() == d_source->indexer()->size() );
}

//---------------------------------------------------------------------------//
// Source physics sends back the number of points it found in its domain
// back to the target.
void Mapper::source_send_point_size()
{
    // Send the number of local points belonging to each target process.
    OrdinateType destination;
    OrdinateType begin_target = 0;
    OrdinateType end_target = d_target->indexer()->size();
    int buffer_size;

    for (int i = begin_target; i < end_target; ++i)
    {
	// Get the global index for the target physics that the buffer is
	// going to.
	destination = d_target->indexer()->l2g(i);

	// Get the number of points in the domain.
	buffer_size = new_map->domain_size(destination);

	// Send a message to the target with the size of the buffer that
	// it will get. 
	nemesis::send_async(&buffer_size, 1, destination);
    }
}

//---------------------------------------------------------------------------//
// Target physics process request for message sizes and post receives.
void Mapper::target_post_receive_buffer(BufferList &buffer_size_list,
					BufferList &buffer_list)
{
    // Initialize.
    OrdinateType src;
    BufferList_Iterator buffer_iter;
    Buffer &buffer;
    int buffer_size;
    denovo::Unpacker u;
	
    // Keep going until all requests have been processed.
    while( !buffer_list.empty() )
    {
	// Find a buffer with a completed communication request.
	buffer_iter = std::find_if(buffer_size_list.begin(), 
				   buffer_size_list.end(), 
				   &Message_Buffer_t::complete);

	// If a completed communication request was found, process it.
	if( buffer_iter != buffer_size_list.end() )
	{
	    // Get the source partition.
	    src = buffer_iter->ordinate();

	    // Get the buffer.
	    buffer = buffer_iter->buffer();

	    // Set the buffer for the unpacker.
	    u.set_buffer( buffer.size(), &buffer[0] );

	    // Get the size of the next buffer we will receive.
	    u >> buffer_size;

	    // Remove this buffer from the size list.
	    buffer_size_list.erase(buffer_iter);

	    // If there is something to receive, post a request.
	    if ( buffer_size > 0 )
	    {
		// Create the buffer and add it to the list.
		buffer_list.push_back( Message_Buffer_t(src, buffer_size) );

		// Clear the buffer just to be safe.
		buffer.clear();

		// Get the request buffer.
		buffer = buffer_list.back();

		// Post asynchronous receive with this receive buffer.
		nemesis::receive_async(buffer.request(),
				       &buffer.buffer()[0],
				       buffer.buffer().size(),
				       buffer.ordinate());
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// Source physics sends its point handles to the targets.
void Mapper::source_send_handles(SP_Transfer_Map new_map)
{
    // For every unique target physics rank in the map, send back the
    // points found in the local domain.
    Set_Iterator destination;
    for (destination = new_map->target_set_begin(); 
	 destination != new_map->target_set_end(); 
	 ++destination)
    {
	// Get the domain iterators for this target rank.
	Map_Pair domain_pair = new_map->domain(*destination);
	Map_Iterator map_it;
	    
	// Create a packer.
	denovo::Packer p;

	// Compute the size of the buffer.
	p.compute_buffer_size_mode();
	for (map_it = domain_pair.first(); 
	     map_it != domain_pair.second();
	     ++map_it)
	{
	    p << (*map_it).second();
	}
	int buffer_size = p.size();

	// Pack the buffer with the handles.
	Buffer buffer(buffer_size);
	p.set_buffer(buffer_size, &buffer[0]);
	for (map_it = domain_pair.first(); 
	     map_it != domain_pair.second();
	     ++map_it)
	{
	    p << (*map_it).second();
	}

	// Send the buffer.
	nemesis::send_async(&buffer[0], buffer_size, *destination);
    }
}

//---------------------------------------------------------------------------//
// Target physics processes handle requests and completes the mapping.
void Mapper::target_process_handles(BufferList &buffer_list,
				    SP_Transfer_Map new_map)
{
    // Initialize.
    OrdinateType src;
    BufferList_Iterator buffer_iter;
    Buffer &buffer;
    int buffer_size;
    int num_handles;
    int j;
    denovo::Unpacker u;
    HandleType handle;

    while ( buffer_list.empty() )
    {
	// Find a buffer with a completed communication request.
	buffer_iter = std::find_if(buffer_list.begin(), 
				   buffer_list.end(), 
				   &Message_Buffer_t::complete);

	// If a completed communication request was found, process it.
	if( buffer_iter != buffer_size_list.end() )
	{
	    // Get the source partition.
	    src = buffer_iter->ordinate();

	    // Get the buffer.
	    buffer = buffer_iter->buffer();

	    // Get the buffer size.
	    buffer_size = buffer.size();

	    // Unpack the buffer.
	    if (buffer_size > 0)
	    {
		// Get the number of handles in the buffer.
		num_handles = buffer_size / sizeof(int);

		// Unpack the handle and put it in the map.
		u.set_buffer(buffer_size, &buffer[0]);

		for (int j = 0; j < num_handles; ++j)
		{
		    u >> handle;

		    new_map->add_target_pair(source, handle);
		}
	    }
	}
    }
}

//---------------------------------------------------------------------------//

// end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Mapper.cc
//---------------------------------------------------------------------------//
