//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   core/Messenger.cc
 * \author Stuart R. Slattery
 * \date   Thu May 26 11:02:57 2011
 * \brief  Messenger class member definitions
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include <algorithm>
#include <vector>

#include "Messenger.hh"
#include "harness/DBC.hh"

namespace coupler
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR and DESTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class DataType_T>
Messenger<DataType_T>::Messenger(const Communicator &comm_global,
				 const std::string &field_name,
				 SP_Physics source,
				 SP_Physics target)
    : d_comm_global(comm_global)
    , d_field_name(field_name)
    , d_source(source)
    , d_target(target)
{
    // Make sure there is a map to operate with.
    Require ( d_source->get_map( d_target->name(), d_field_name ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class DataType_T>
Messenger<DataType_T>::~Messenger()
{ /* ... */ }

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Communicate the field from the source to the target.
 */
template<class DataType_T>
void Messenger<DataType_T>::communicate()
{
    // Set the internal communicator.
    nemesis::set_internal_comm(d_comm_global);

    // Create an empty list of message buffers.
    BufferList buffer_list;

    // Target physics posts receives.
    if ( d_target->te() )
    {
	post_receives(buffer_list);
    }

    // Source physics sends buffers to target physics.
    if ( d_source->te() )
    {
	send();
    }

    // Target physics processes buffers from source physics.
    if ( d_target->te() )
    {
	process_requests(buffer_list);
    }

    // Reset the internal communicator.
    nemesis::reset_internal_comm();
}


//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Post asynchronous receives for the buffers.
 * \param buffer_list List of buffers containing data.
 */
template<class DataType_T>
void Messenger<DataType_T>::post_receives(BufferList &buffer_list)
{
    // Initialize.
    Message_Buffer_t &buffer;
    int buffer_size = 0;

    // Create the buffers and post the receives for each source process
    // sending to this target process.
    Set_Iterator src;
    Set_Pair src_bound = 
	d_source->get_map( d_target->name(), d_field_name )->sources(); 

    for ( src = src_bound.first(); src != src_bound.second(); ++src) 
    {
	Check ( *src < nemesis::nodes() );

        // Compute the size of the buffer.
	buffer_size = d_source->get_map( 
	    d_target->name(), d_field_name )->range_size(src) 
		      * ( sizeof(HandleType) + sizeof(DataType) );

	// Create the buffer and add it to the list.
        buffer_list.push_back( Message_Buffer_t(*src, buffer_size) );

        // Get the request buffer.
        buffer = buffer_list.back();

        // Post asynchronous receive with this receive buffer.
        nemesis::receive_async(buffer.request(), 
			       &buffer.buffer()[0], 
                               buffer.buffer().size(), 
			       buffer.ordinate());
    }
    
    // Make sure we made all of the buffers we're going to receive.
    Ensure ( buffer_list.size() == 
	     std::distance( src_bound.first(), src_bound.second() ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do asynchronous sends of data from source to target.
 */
template<class DataType_T>
void Messenger<DataType_T>::send()
{
    // Initialize.
    denovo::Packer p;

    std::vector<DataType> data;
    typename std::vector<DataType>::const_iterator data_it;

    std::vector<HandleType> handles;
    typename std::vector<HandleType>::const_iterator handle_it;

    Buffer buffer;
    int buffer_size;

    Map_Iterator domain_it;
    Map_Pair domain_bound;

    // Loop over the target partitions and send the data.    
    Set_Iterator destination;
    Set_Pair destination_bound = 
	d_source->get_map( d_target->name(), d_field_name )->targets();

    for ( destination = destination_bound.first();
	  destination != destination_bound.second();
	  ++destination)     
    {
        Check ( *destination < nemesis::nodes() );
	
	// Clear the buffer.
	buffer.clear();

        // Resize the buffer.
	buffer_size = d_source->get_map( 
	    d_target->name(), d_field_name )->domain_size(*destination) 
		      * ( sizeof(HandleType) + sizeof(DataType) );

	buffer.resize(buffer_size);

        // Register the current buffer with the packer
        p.set_buffer( buffer.size(), &buffer[0] );

        // Pack the data we pull from the source physics.
	domain_bound = d_source->get_map( 
	    d_target->name(), d_field_name )->domain(*destination);

	for (domain_it = domain_bound.first(); 
	     domain_it != domain_bound.second();
	     ++domain_it)
        {
	    handles.push_back( domain_it->second() );
	}

	data = d_source->te()->pull_data(d_field_name, handles, data);
	Check ( handles.size() == data.size() );

	for (data_it = data.begin(), handle_it = handles.begin();
	     data_it != data.end();
	     ++data_it, ++handle_it)
        {
	    p << *handle_it;
	    p << *data_it;
	}

        // blocking send the buffer to the remote process
        nemesis::send_async( &buffer[0], buffer.size(), *destination );

	// Clear the vectors.
	handles.clear();
	data.clear();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process the target requests and push the data onto the targets.
 * \param buffer_list List of buffers containing data.
 */
template<class DataType_T>
void Messenger<DataType_T>::process_requests(BufferList &buffer_list)
{
    // Initialize.
    denovo::Unpacker u;

    HandleType temp_handle;
    std::vector<HandleType> handles;
    typename std::vector<HandleType>::const_iterator handle_it;

    DataType temp_data;
    std::vector<DataType> data;
    typename std::vector<DataType>::const_iterator data_it;

    OrdinateType src;

    BufferList_Iterator buffer_iter;
    Buffer &buffer;

    Map_Iterator range_it;
    Map_Pair range_bound;
    
    // Keep going until all requests have been processed.
    while( !buffer_list.empty() )
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

            // Set the buffer for the unpacker.
            u.set_buffer( buffer.size(), &buffer[0] );
            
            // Loop over the target handles and push the data.
	    range_bound = d_source->get_map( 
		d_target->name(), d_field_name )->range(src);

            for( range_it = range_bound.first();
		 range_it != range_bound.second();
		 ++range_it)
            {
		u >> temp_handle;
                u >> temp_data;
		
		handles.push_back( temp_handle );
		data.push_back( temp_data );
	    }
	    
	    Check ( handles.size() == data.size() );
	    d_target->te()->push_data(d_field_name, handles, data);

            // Remove this buffer from the list.
            buffer_list.erase(buffer_iter);

	    // Clear the vectors.
	    handles.clear();
	    data.clear();
        }
    }
}

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Messenger.cc
//---------------------------------------------------------------------------//
