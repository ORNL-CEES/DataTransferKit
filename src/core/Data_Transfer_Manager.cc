//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Data_Transfer_Manager.cc
 * \author Stuart Slattery
 * \date   Wed Oct 05 11:02:44 2011
 * \brief  Data_Transfer_Manager member definitons.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "Data_Transfer_Manager.hh"
#include "harness/DBC.hh"
#include <vector>
#include <iterator>

namespace coupler
{

//---------------------------------------------------------------------------//
// Constructor.
Data_Transfer_Manager::Data_Transfer_Manager(const Communicator &comm_global)
    : d_comm_global(comm_global)
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
Data_Transfer_Manager::~Data_Transfer_Manager()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Register a physics to be controlled by the manager.
void Data_Transfer_Manager::add_physics(const std::string &physics_name,
					Transfer_Evaluator *te)
{
    // Make a physics object.
    SP_Physics new_physics = new Physics(te, d_comm_global);

    // Add it to the physics database.
    d_physics_db.insert( 
	std::pair<std::string,SP_Physics>(physics_name, new_physics) );
}

//---------------------------------------------------------------------------//
// Register a field to be controlled by the manager.
void Data_Transfer_Manager::add_field(const std::string &field_name)
{
    // Require that the field is supported by the each code.
    Require( d_te_a->register_field(field_name) );
    Require( d_te_b->register_field(field_name) );
}

//---------------------------------------------------------------------------//
// Build the topology map for transfer from a source physics to a target
// physics for a particular field.
void Data_Transfer_Manager::map(const std::string &field_name,
				const std::string &source_physics,
				const std::string &target_physics);
{
    // Operate on the global communicator.
    nemesis::set_internal_comm(d_comm_global);

    // Make a transfer map.
    SP_Transfer_Map new_map = new Transfer_Map();

    // Get the transfer evaluator implementations for the physics we are
    // operating on.
    SP_Physics source = d_physics_db[source_physics];
    SP_Physics target = d_physics_db[target_physics];

    // Set the iteration bounds for loops over the source and target process
    // ids.
    int begin_source = 0;
    int end_source = 0;
    int begin_target = 0;
    int end_target = 0;
    if( source->te() || target->te() )
    {
	end_source = source->indexer()->size();
	end_target = target->indexer()->size();
    }

    // Target point coordinate vector iterators.
    Coord_Iterator points_begin, points_end;

    // Target point handle vector iterators.
    Handle_Iterator handles_begin, handles_end;

    // Target physics registers its target points for the field being mapped.
    target->te()->register_xyz(field_name, 
			       points_begin, points_end, 
			       handles_begin, handles_end);
    Check( std::distance(points_begin, points_end) % 3 == 0 );
    Check( std::distance(points_begin, points_end) / 3 == 
	   std::distance(handles_begin, handles_end) );

    // Target physics sends all of its target points to each source physics
    // process.
    if (target->te())
    {
	// Build a buffer of the local points to send to the source physics.
	Buffer buffer;

	// Create a packer.
	denovo::Packer p;

	// Compute the size of the buffer.
	Coord_Iterator coord_iter = points_begin;
	Handle_Iterator handle_iter;
	p.compute_buffer_size_mode();
	for (handle_iter = handles_begin, handle_iter != handels_end; ++iter)
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
	}
	int buffer_size = buffer.size();

	// Send the local target points to all processes of the source physics.
	int destination;
	for (int i = begin_source; i < end_source; ++i)
	{
	    // Get the global index for the source physics that the buffer is
	    // going to.
	    destination = source->indexer()->l2g(i);

	    // Send a message to the source with the size of the buffer that it will
	    // get. 
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
    if (source->te())
    {
	// Source will get a message with target points from every target process.
	for (int i = begin_target; i < end_target; ++i)
	{
	    // Get the global index of the target process.
	    int src = target->indexer()->l2g(i);

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
		Handle handle;
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
		    if ( source->te()->find_xyz(x, y, z, handle) )
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
    if (source->te())
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
    if (target->te())
    {
	// Get a buffer from all source processes.
	for (int i = begin_source; i < end_source; ++i)
	{
	    // Get the global rank of A we are receiving from.
	    int src = source->indexer()->l2g(i);

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
		    Handle handle;
		    u >> handle;

		    new_map->add_target_pair(source, handle);
		}
	    }
	}
    }

    // Add the new map to the source physics database.
    source->set_map(target_physics, field_name, new_map);
}

//---------------------------------------------------------------------------//
// Transfer data associated with a field from a source physics to a target
// physics. 
void Data_Transfer_Manager::transfer(const std::string &field_name,
				     const std::string &source_physics,
				     const std::string &target_physics)
{

}

//---------------------------------------------------------------------------//

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Data_Transfer_Manager.cc
//---------------------------------------------------------------------------//
