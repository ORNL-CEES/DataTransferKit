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
Messenger<DataType_T>::Messenger()
{
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
void Messenger<DataType_T>::communicate(
    const Communicator &comm_global,
    SP_Transfer_Data_Field transfer_data_field)
{
    // Set the internal communicator.
    nemesis::set_internal_comm(comm_global);

    // Create an empty list of message buffers.
    BufferList buffer_list;

    // Target physics posts receives.
    if ( transfer_data_field->target() )
    {
	post_receives(buffer_list);
    }

    // Source physics sends buffers to target physics.
    if ( transfer_data_field->source() )
    {
	send();
    }

    // Target physics processes buffers from source physics.
    if ( transfer_data_field->target() )
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
void Messenger<DataType_T>::post_receives(
    SP_Transfer_Data_Field transfer_data_field,
    BufferList &buffer_list)
{
    // Initialize.
    Message_Buffer_t &buffer;
    int buffer_size = 0;

    // Create the buffers and post the receives for each source process
    // sending to this target process.
    Set_Iterator src;
    Set_Pair src_bound = transfer_data_field->get_map()->sources(); 

    for ( src = src_bound.first(); src != src_bound.second(); ++src) 
    {
	Check ( *src < nemesis::nodes() );

        // Compute the size of the buffer.
	buffer_size = transfer_data_field->get_map()->range_size(src)
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
void Messenger<DataType_T>::send(SP_Transfer_Data_Field transfer_data_field)
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
    Set_Pair destination_bound = transfer_data_field->get_map()->targets();

    for ( destination = destination_bound.first();
	  destination != destination_bound.second();
	  ++destination)     
    {
        Check ( *destination < nemesis::nodes() );
	
	// Clear the buffer.
	buffer.clear();

        // Resize the buffer.
	buffer_size = transfer_data_field->get_map()->domain_size(*destination)
		      * ( sizeof(HandleType) + sizeof(DataType) );

	buffer.resize(buffer_size);

        // Register the current buffer with the packer
        p.set_buffer( buffer.size(), &buffer[0] );

        // Get the handles we will send.
	domain_bound = transfer_data_field->get_map()->domain(*destination);

	for (domain_it = domain_bound.first(); 
	     domain_it != domain_bound.second();
	     ++domain_it)
        {
	    handles.push_back( domain_it->second() );
	}

	// Get the data we will send.
	data = transfer_data_field->source()->send_data( 
	    transfer_data_field->name(), handles, data);
	Check ( handles.size() == data.size() );

	// Pack the handles and data into a buffer.
	for (data_it = data.begin(), handle_it = handles.begin();
	     data_it != data.end();
	     ++data_it, ++handle_it)
        {
	    p << *handle_it;
	    p << *data_it;
	}

        // non-blocking send the buffer to the remote process
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
void Messenger<DataType_T>::process_requests(
    SP_Transfer_Data_Field transfer_data_field,
    BufferList &buffer_list)
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
            
            // Unpack the data from the buffer.
	    range_bound = transfer_data_field->get_map()->range(src);

            for( range_it = range_bound.first();
		 range_it != range_bound.second();
		 ++range_it)
            {
		u >> temp_handle;
                u >> temp_data;
		
		handles.push_back( temp_handle );
		data.push_back( temp_data );
	    }
	    
	    // Push the data onto the targets.
	    Check ( handles.size() == data.size() );
	    transfer_data_field->target()->receive_data(
		transfer_data_field->name(), handles, data);

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
