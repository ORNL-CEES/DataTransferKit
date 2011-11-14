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
 * \param comm_global The global communicator encapsulating the all processes
 * participating in coupling.
 * \param field_name The name of the field being mapped.
 * \param source Smart pointer to source physics.
 * \param target Smart pointer to target physics.
 */
template<class DataType_T>
Mapper<DataType_T>::Mapper(const Communicator &comm_global,
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
template<class DataType_T>
Mapper<DataType_T>::~Mapper()
{ /* ... */ }


//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Map the field from the source onto the target.
 */
template<class DataType_T>
void Mapper<DataType_T>::map()
{
    //  Set the internal communicator.
    nemesis::set_internal_comm(d_comm_global);

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
	source_send_point_size(new_map);
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
    d_source->set_map(d_target->name(), d_field_name, new_map);
    
    // Reset the internal communicator.
    nemesis::reset_internal_comm();
}


//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Source physics post receives for buffer sizes.
 * \param buffer_size_list A buffer list of buffers containing the number of
 * points that will come from each target process to this source process.
 */
template<class DataType_T>
void Mapper<DataType_T>::source_post_receive_size(
    BufferList &buffer_size_list)
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
/*!
 * \brief Target physics sends point sizes to source.
 * \param points_begin Iterator to the beginning of the point coordinate
 * array.
 * \param points_end Iterator to the end of the point coordinate array. 
 * \param handles_begin Iterator to the beginning of the point handle array.
 * \param handles_end Iterator to the end of the point handle array.
 */
template<class DataType_T>
void Mapper<DataType_T>::target_send_point_size(Coord_Iterator &points_begin,
						Coord_Iterator &points_end,
						Handle_Iterator &handles_begin,
						Handle_Iterator &handles_end)
{
    // Target physics registers its target points for the field being mapped.
    d_target->te()->register_xyz(d_field_name, 
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
/*!
 * \brief Source physics process requests for message sizes and post receives
 * for buffers. 
 * \param buffer_size_list Buffer list with buffers holding the number of
 * points from each target process.
 * \param buffer_list Buffer list with buffers holding points from each target
 * process.
 */
template<class DataType_T>
void Mapper<DataType_T>::source_post_receive_buffer(
    BufferList &buffer_size_list,
    BufferList &buffer_list)
{
    // Initialize.
    OrdinateType src;
    BufferList_Iterator buffer_iter;
    Message_Buffer_t &buffer;
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
	    buffer = buffer_iter;

	    // Set the buffer for the unpacker.
	    u.set_buffer( buffer.buffer().size(), &buffer.buffer()[0] );

	    // Get the size of the next buffer we will receive.
	    u >> buffer_size;

	    // Remove this buffer from the size list.
	    buffer_size_list.erase(buffer_iter);

	    // Create the buffer and add it to the list.
	    buffer_list.push_back( Message_Buffer_t(src, buffer_size) );

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
/*!
 * \brief Target send points to source.
 * \param points_begin Iterator to the beginning of the point coordinate
 * array.
 * \param points_end Iterator to the end of the point coordinate array. 
 * \param handles_begin Iterator to the beginning of the point handle array.
 * \param handles_end Iterator to the end of the point handle array.
 */
template<class DataType_T>
void Mapper<DataType_T>::target_send_points(Coord_Iterator points_begin,
					    Coord_Iterator points_end,
					    Handle_Iterator handles_begin,
					    Handle_Iterator handles_end)
{
    // Build a buffer of the local points to send to the source physics.
    Buffer buffer;
    int buffer_size;
    denovo::Packer p;
    
    // Compute the size of the buffer.
    Coord_Iterator coord_iter = points_begin;
    Handle_Iterator handle_iter;
    p.compute_buffer_size_mode();
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
    buffer_size = p.size();

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
    buffer_size = buffer.size();

    // Send the local target points to all processes of the source physics.
    OrdinateType destination;
    OrdinateType begin_source = 0;
    OrdinateType end_source = d_source->indexer()->size();

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
/*!
 * \brief Source physics process request and build part of the map.
 * \param buffer_list Buffer list of buffers containing points from target
 * processes. 
 * \param new_map Smart pointer to the transfer map being generated by the
 * mapping algorithm.
 */
template<class DataType_T>
void Mapper<DataType_T>::source_process_points(BufferList &buffer_list,
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
    CoordinateType x, y, z;

    while ( buffer_list.empty() )
    {
	// Find a buffer with a completed communication request.
	buffer_iter = std::find_if(buffer_list.begin(), 
				   buffer_list.end(), 
				   &Message_Buffer_t::complete);

	// If a completed communication request was found, process it.
	if( buffer_iter != buffer_list.end() )
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
/*! 
 * \brief Target physics post receives for return buffer size.
 * \param buffer_size_list Buffer list of buffers containing the number of
 * handles from each source process.
 */
template<class DataType_T>
void Mapper<DataType_T>::target_post_receive_size(BufferList &buffer_size_list)
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
/*!
 * \brief Source physics sends back the number of points it found in its
 * domain back to the target.
 */
template<class DataType_T>
void Mapper<DataType_T>::source_send_point_size(SP_Transfer_Map new_map)
{
    // Send the number of local points belonging to each target process.
    int buffer_size;
    OrdinateType destination;
    OrdinateType begin_target = 0;
    OrdinateType end_target = d_target->indexer()->size();

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
/*!
 * \brief Target physics process request for message sizes and post receives.
 * \param buffer_size_list Buffer list with buffers holding the number of
 * points from each target process.
 * \param buffer_list Buffer list with buffers holding points from each target
 * process.
 */
template<class DataType_T>
void Mapper<DataType_T>::target_post_receive_buffer(
    BufferList &buffer_size_list,
    BufferList &buffer_list)
{
    // Initialize.
    OrdinateType src;
    BufferList_Iterator buffer_iter;
    Message_Buffer_t &buffer;
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
	    buffer = buffer_iter;

	    // Set the buffer for the unpacker.
	    u.set_buffer( buffer.buffer().size(), &buffer.buffer()[0] );

	    // Get the size of the next buffer we will receive.
	    u >> buffer_size;

	    // Remove this buffer from the size list.
	    buffer_size_list.erase(buffer_iter);

	    // If there is something to receive, post a request.
	    if ( buffer_size > 0 )
	    {
		// Create the buffer and add it to the list.
		buffer_list.push_back( Message_Buffer_t(src, buffer_size) );

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
/*! 
 * \brief Source physics sends its point handles to the targets.
 * \param new_map Smart pointer to the Transfer_Map being generated by the
 * mapping algorithm.
 */
template<class DataType_T>
void Mapper<DataType_T>::source_send_handles(SP_Transfer_Map new_map)
{
    // For every unique target physics rank in the map, send back the
    // points found in the local domain.
    int buffer_size;
    denovo::Packer p;
    Set_Iterator destination;
    Set_Pair destination_bound 	= 
	d_source->get_map( d_target->name(), d_field_name )->targets();

    for (destination = destination_bound.first(); 
	 destination != destination_bound.second(); 
	 ++destination)
    {
	// Get the domain iterators for this target rank.
	Map_Pair domain_pair = new_map->domain(*destination);
	Map_Iterator map_it;
	    
	// Compute the size of the buffer.
	p.compute_buffer_size_mode();
	for (map_it = domain_pair.first(); 
	     map_it != domain_pair.second();
	     ++map_it)
	{
	    p << (*map_it).second();
	}
	buffer_size = p.size();

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
/*!
 * \brief Target physics processes handle requests and completes the mapping.
 * \param buffer_list Buffer list of buffers containing points from target
 * processes. 
 * \param new_map Smart pointer to the transfer map being generated by the
 * mapping algorithm.
 */
template<class DataType_T>
void Mapper<DataType_T>::target_process_handles(BufferList &buffer_list,
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
	if( buffer_iter != buffer_list.end() )
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

		    new_map->add_range_pair(src, handle);
		}
	    }
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Mapper.cc
//---------------------------------------------------------------------------//
