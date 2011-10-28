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

namespace coupler
{

//---------------------------------------------------------------------------//
// Constructor.
Data_Transfer_Manager::Data_Transfer_Manager(Communicator comm_global,
					     Transfer_Evaluator* TE_A,
					     Transfer_Evaluator* TE_B)
    : d_comm_global(comm_global)
{
    // Wrap the raw Transfer_Evaluator pointers.
    d_te_a = TE_A;
    d_te_b = TE_B;

    // Get the physics' subcommunicators.
    d_te_a->register_comm(d_comm_a);
    d_te_b->register_comm(d_comm_b);

    // Generate local to global rank indexers.
    if(d_te_a || d_te_b)
    {
	d_indexer_a = new LG_Indexer(d_comm_global, d_comm_a, d_te_a);
	Ensure( d_indexer_a );

	d_indexer_b = new LG_Indexer(d_comm_global, d_comm_b, d_te_b);
	Ensure( d_indexer_b );
    }
}

//---------------------------------------------------------------------------//
// Destructor.
Data_Transfer_Manager::~Data_Transfer_Manager()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Register a field to be controlled by the manager.
void Data_Transfer_Manager::add_field(std::string field_name)
{
    // Check that the field is supported by the each code.
    Require( d_te_a->register_field(field_name) );
    Require( d_te_b->register_field(field_name) );

    // Add the new field to the database.
    d_f_db.add_field(field_name);
}

//---------------------------------------------------------------------------//
// Build the topology map for transfer from A to B for a particular field. B
// provides the target points and A is the source.
void Data_Transfer_Manager::map_A2B(std::string field_name)
{
    // Operate on the global communicator.
    nemesis::set_internal_comm(d_comm_global);

    // Initialize map.
    d_map_a2b = new Transfer_Map();    

    // Set the iteration bounds for loops over the A and B process ids.
    int begin_a = 0;
    int end_a = 0;
    int begin_b = 0;
    int end_b = 0;
    if( d_te_a || d_te_b )
    {
	end_a = d_indexer_a->size();
	end_b = d_indexer_b->size();
    }

    // Target point coordinate vector iterators.
    Coord_Iterator points_begin, points_end;

    // Target point handle vector iterators.
    Handle_Iterator handles_begin, handles_end;

    // B registers its target points.
    d_te_b->register_xyz(field_name, 
			 points_begin, points_end, 
			 handles_begin, handles_end);
    Check( std::distance(points_end, points_begin) % 3 == 0 );
    Check( std::distance(points_end,points_begin) / 3 == 
	   std::distance(handles_begin, handles_end) );

    // B sends all of its target points to each A process.
    if (d_te_b)
    {
	// Build a buffer of the local points to send to A.
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
	}

	int buffer_size = buffer.size();

	// Send the local points to all processes of A.
	int destination;
	for (int i = begin_a; i < end_a; ++i)
	{
	    // Get the global index for A that the buffer is going to.
	    destination = d_indexer_a->l2g(i);

	    // Send a message to A with the size of the buffer that it will
	    // get. 
	    nemesis::send_async(&buffer_size, 1, destination);

	    // Send the buffer of points.
	    if (buffer_size > 0)
	    {
		nemesis::send_async(&buffer[0], buffer_size, destination);
	    }
	}
    }

    // A receives all of the target points from A and builds the topology
    // map. 
    if (d_te_a)
    {
	// A will get a message with target points from every B process.
	for (int i = begin_b; i < end_b; ++i)
	{
	    // Get the global index of the B process.
	    int source = d_indexer_b->l2g(i);

	    // Get the size of the incoming target point buffer.
	    int buffer_size;
	    nemesis::receive(&buffer_size, 1, source);

	    // Unpack the points and add them to the map.
	    int num_points;
	    if (buffer_size > 0)
	    {
		// Receive the buffer from B.
		Buffer buffer(buffer_size);
		nemesis::receive(&buffer[0], buffer_size, source);

		// Compute the number of points in the buffer.
		num_points = buffer_size / ( sizeof(double) * 3 );

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

		    // See if this point is in the spatial domain of A.
		    if ( d_te_a->find_xyz(x, y, z, handle) )
		    {
			// Add the handle to the map with the target rank.
			d_map_a2b->add_domain_pair(source, handle);
		    }
		}
	    }
	}
    }

    // Barrier after sending target points from B to A.
    nemesis::global_barrier();

    // Send all target points found in A back to B to complete the map.
    if (d_te_a)
    {
	// For every unique B in the map, send back the points found in the
	// local domain.
	Set_Iterator destination;
	for (destination = d_map_a2b->target_set_begin(); 
	     destination != d_map_a2b->target_set_end(); 
	     ++destination)
	{
	    // Get the domain iterators for this target rank.
	    Iterator_Pair domain_pair = map_a2b->domain(*destination);
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

    // B gets a message from A with the handles to its target points that A
    // found. B associates these handles with the rank of A that they came
    // from to complete the map.
    if (d_te_b)
    {
	// Get a buffer from all A processes.
	for (int i = begin_a; i < end_a; ++i)
	{
	    // Get the global rank of A we are receiving from.
	    int source = d_indexer_a->l2g(i);

	    // Receive the buffer size.
	    int buffer_size;
	    nemesis::receive(&bufer_size, 1, source);

	}
    }
}

//---------------------------------------------------------------------------//
// Transfer data from A to B.
void Data_Transfer_Manager::transfer_A2B(std::string field_name)
{

}

//---------------------------------------------------------------------------//

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Data_Transfer_Manager.cc
//---------------------------------------------------------------------------//
