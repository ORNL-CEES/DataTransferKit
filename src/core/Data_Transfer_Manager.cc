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
#include <algorithm>

namespace coupler
{

//---------------------------------------------------------------------------//
// Constructor.
Data_Transfer_Manager::Data_Transfer_Manager(Communicator_t comm_global,
					     Transfer_Evaluator* TE_A,
					     Transfer_Evaluator* TE_B)
    : d_comm_global(comm_global)
{
    // operate on the global communicator
    nemesis::set_internal_comm(d_comm_global);

    // Wrap the raw pointers.
    d_te_a = TE_A;
    d_te_b = TE_B;

    // Get the physics' subcommunicators.
    d_te_a->register_comm(d_comm_a);
    d_te_b->register_comm(d_comm_b);

    // Generate local to global indexers.
    if(d_te_a || d_te_b)
    {
	d_indexer_a = new LG_Indexer(d_comm_global, d_comm_a, d_te_a);
	Ensure( d_indexer_a );

	d_indexer_b = new LG_Indexer(d_comm_global, d_comm_b, d_te_b);
	Ensure( d_indexer_b );
    }

    // reset the internal communicator
    nemesis::reset_internal_comm();
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
    // operator on the global communicator
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

    // Target point coordinate vector.
    std::vector<double> points;

    // Target point handle vector.
    std::vector<Handle> handles;

    // B registers its target points.
    d_te_b->register_xyz(field_name, points, handles);
    Check( points.size() % 3 == 0 );
    Check( points.size() / 3 == handles.size() );

    // B sends all of its target points to each A process.
    if (d_te_b)
    {
	// Build a buffer of the local points to send to A.
	Buffer buffer;
	build_buffer(buffer, points.begin(), points.end());
	int buffer_size = buffer.size();

	// Send the local points to all processes of A.
	for (int i = begin_a; i < end_a; ++i)
	{
	    // Get the global index for A that the buffer is going to.
	    int destination = d_indexer_a->l2g(i);

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
	    if (buffer_size > 0)
	    {
		// Receive the buffer from B.
		Buffer buffer(buffer_size);
		nemesis::receive(&buffer[0], buffer_size, source);

		// Compute the number of points in the buffer.
		int num_points = buffer_size / 
				 ( sizeof(double) * 3 );

		// Unpack the buffer.
		denovo::Unpacker u;
		u.set_buffer(buffer_size, &buffer[0]);
		for (int j = 0; j < num_points; ++j)
		{
		    double x, y, z;
		    u >> x;
		    u >> y;
		    u >> z;

		    // See if this point is in the spatial domain of A.
		    Handle local_handle;
		    if ( d_te_a->find_xyz(x, y, z, local_handle) )
		    {
			// Add the local handle to the map with the target
			// rank.
			map_a2b->add_domain_pair(source, local_handle);
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
	
    }
}

//---------------------------------------------------------------------------//
// Transfer data from A to B.
void Data_Transfer_Manager::transfer_A2B(std::string field_name)
{
    // Get the field we are operating on.
    d_f_db;

    // Get the map contents.
    Transfer_Map::Vector_Int index = d_map_A2B->get_index();
    Transfer_Map::Vector_Int rank = d_map_A2B->get_rank();

    // Register the domain of A. We want to do this every time because it may
    // have changed.
    Data_Iterator domain_begin;
    Data_Iterator domain_end;
    d_te_a->register_domain(field_name, domain_begin, domain_end);
    int domain_size = domain_end - domain_begin;
    
    // Register the range of B. We want to do this every time because it may
    // have changed.
    Iterator range_begin;
    Iterator range_end;
    if (d_te_a->register_range(field_name, range_begin, range_end) )
    int range_size = range_end - range_begin;

    // Transfer A to B.
    std::vector<ValueType> data(domain_size);

    // Pull the data from A.
    std::copy(domain_begin, domain_end, data.begin());
    
    // Modify it with the map before pushing it to B.
    

    // Check the data vector size before copy.
    assert( data.size() == range_end - range_begin );
    
    // Push the data to B.
    std::copy(data.begin(), data.end(), range_begin);
}

//---------------------------------------------------------------------------//
// Generate a buffer and pack it with data.
template<class X>
void Data_Transfer_Manager::build_buffer(Buffer &buffer, 
					 X::value_type::const_iterator begin,
					 X::value_type::const_iterator end)
{
    // Clear the buffer to make sure it's empty.
    buffer.clear();

    // Create a packer.
    denovo::Packer p;

    // Data iterator.
    typename X::value_type::const_iterator iter;

    // Compute the size of the buffer.
    p.compute_buffer_size_mode();
    for (iter = begin, iter != end; ++iter)
    {
	p << *iter;
    }
    int buffer_size = p.size();

    // Set the size of the buffer.
    buffer.resize(buffer_size);

    // Pack the data into the buffer.
    Check( buffer.size() == buffer_size );
    if (buffer_size > 0)
    {
	p.set_buffer(buffer_size, &buffer[0]);
	for (iter = begin, iter != end; ++iter)
	{
	    p << *iter;
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Data_Transfer_Manager.cc
//---------------------------------------------------------------------------//
