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
    BufferList buffer_list;

    // Source physics post receives.
    if ( d_source->te() )
    {
	source_post_receives(buffer_list);
    }

    // Target physics sends to the source.
    if ( d_target->te() )
    {
	target_send();
    }

    // Source physics processes buffers from target and adds to the map.
    if ( d_source->te() )
    {
	source_process_requests(buffer_list, new_map);
    }

    // Barrier before continuing.
    nemesis::global_barrier();
    buffer_list.clear();

    // Target physics post receives.
    if ( d_target->te() )
    {
	target_post_receives(buffer_list);
    }

    // Source physics sends to the target.
    if ( d_source->te() )
    {
	source_send();
    }

    // Target physics processes buffers from source and finishes the map.
    if ( d_target->te() )
    {
	target_process_requests(buffer_list, new_map);
    }
    
    // Reset the internal communicator.
    nemesis::reset_internal_comm();
}




    // Set the iteration bounds for loops over the source and target process
    // ids.
    int begin_source = 0;
    int end_source = 0;
    int begin_target = 0;
    int end_target = 0;
    if( d_source->te() || d_target->te() )
    {
	end_source = d_source->indexer()->size();
	end_target = d_target->indexer()->size();
    }


    // Target physics sends all of its target points to each source physics
    // process.
    if (d_target->te())
    {
	// Target point coordinate vector iterators.
	Coord_Iterator points_begin, points_end;

	// Target point handle vector iterators.
	Handle_Iterator handles_begin, handles_end;

	// Target physics registers its target points for the field being mapped.
	d_target->te()->register_xyz(field_name, 
				     points_begin, points_end, 
				     handles_begin, handles_end);

	Check( std::distance(points_begin, points_end) % 3 == 0 );
	Check( std::distance(points_begin, points_end) / 3 == 
	       std::distance(handles_begin, handles_end) );

	// Build a buffer of the local points to send to the source physics.
	Buffer buffer;

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

	    // Send a message to the source with the size of the buffer that
	    // it will get. 
	    nemesis::send_async(&buffer_size, 1, destination);

	    // Send the buffer of points.
	    if (buffer_size > 0)
	    {
		nemesis::send_async(&buffer[0], buffer_size, destination);
	    }
	}
    }

    // Source physics receives all of the target points from the target
    // physics and builds the topology map. 
    if (d_source->te())
    {
	// Source will get a message with target points from every target process.
	for (int i = begin_target; i < end_target; ++i)
	{
	    // Get the global index of the target process.
	    int src = d_target->indexer()->l2g(i);

	    // Get the size of the incoming target point buffer.
	    int buffer_size;
	    nemesis::receive(&buffer_size, 1, src);

	    // Unpack the points and add them to the map.
	    int num_points;
	    if (buffer_size > 0)
	    {
		// Receive the buffer from the target process.
		Buffer buffer(buffer_size);
		nemesis::receive(&buffer[0], buffer_size, src);

		// Compute the number of points in the buffer.
		num_points = buffer_size / ( sizeof(double) * 3 + sizeof(int));

		// Unpack the buffer.
		HandleType handle;
		Coordinate x, y, z;
		denovo::Unpacker u;
		u.set_buffer(buffer_size, &buffer[0]);
		for (int j = 0; j < num_points; ++j)
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
	    }
	}
    }

    // Barrier after sending target points from the target physics to the
    // source physics.
    nemesis::global_barrier();

    // Send all target points found in the source physics back to the target
    // physics to complete the map.
    if (d_source->te())
    {
	// For every unique target physics rank in the map, send back the
	// points found in the local domain.
	Set_Iterator destination;
	for (destination = new_map->target_set_begin(); 
	     destination != new_map->target_set_end(); 
	     ++destination)
	{
	    // Get the domain iterators for this target rank.
	    Iterator_Pair domain_pair = new_map->domain(*destination);
	    std::multimap<int,int>::const_iterator mapit;
	    
	    // Create a packer.
	    denovo::Packer p;

	    // Compute the size of the buffer.
	    p.compute_buffer_size_mode();
	    for (map_it = domain_pair.first(); 
		 map_it != domain_pair.second;
		 ++map_it)
	    {
		p << (*map_it).second;
	    }
	    int buffer_size = p.size();

	    // Pack the buffer.
	    Buffer buffer(buffer_size);
	    p.set_buffer(buffer_size, &buffer[0]);
	    for (map_it = domain_pair.first(); 
		 map_it != domain_pair.second;
		 ++map_it)
	    {
		p << (*map_it).second;
	    }

	    // Send the size of the buffer.
	    nemesis::send_async(&buffer_size, 1, *destination);

	    // Send the buffer.
	    nemesis::send_async(&buffer[0], buffer_size, *destination);
	}
    }

    // The target physics gets a message from the source physics with the
    // handles to its target points that the source found. The target
    // associates these handles with the rank of the source that they came
    // from to complete the map. 
    if (d_target->te())
    {
	// Get a buffer from all source processes.
	for (int i = begin_source; i < end_source; ++i)
	{
	    // Get the global rank of A we are receiving from.
	    int src = d_source->indexer()->l2g(i);

	    // Receive the buffer size.
	    int buffer_size;
	    nemesis::receive(&buffer_size, 1, src);

	    // Unpack the buffer.
	    if (buffer_size > 0)
	    {
		// Create a buffer.
		Buffer buffer(buffer_size);

		// Receive the buffer.
		nemesis::receive(&buffer[0], buffer_size, src);
		
		// Get the number of handles in the buffer.
		int num_handles = buffer_size / sizeof(int);

		// Unpack the handle and put it in the map.
		denovo::Unpacker u;
		u.set_buffer(buffer_size, &buffer[0]);
		for (int j = 0; j < num_handles; ++j)
		{
		    HandleType handle;
		    u >> handle;

		    new_map->add_target_pair(source, handle);
		}
	    }
	}
    }

    // Add the new map to the source physics database.
    d_source->set_map(target_physics, field_name, new_map);
} 


//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
// Source physics post receives.
void source_post_receives()
{

}

// Target physics send to source.
void target_send()
{

}

// Source physics process requests.
void source_process_requests()
{

}
    
// Target physics post receives.
void target_post_receives()
{

}

// Source physics send to target.
void source_send()
{

}

// Target physics process requests.
void source_process_requests()
{

}


// end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Mapper.cc
//---------------------------------------------------------------------------//
