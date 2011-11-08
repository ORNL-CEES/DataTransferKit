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

#include "Messenger.hh"
#include "harness/DBC.hh"

namespace coupler
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Messenger::Messenger(const Communicator &comm_global,
		     const std::string &field_name,
		     SP_Physics source,
		     SP_Physics target)
    : d_comm_global(comm_global)
    , d_field_name(field_name)
    , d_source(source)
    , d_target(target)
{  
    // Make sure there is a map to operate with.
    Ensure ( d_source->get_map( d_target->name(), field_name ) );

    // Create the buffer vector and compute their sizes
    calculate_buffer_sizes();
}


//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Transfer data.
 */
void Messenger::communicate()
{
    // Set the internal communicator.
    nemesis::set_internal_comm(d_comm_global);

    // Create an empty list of message buffers.
    std::list<Message_Buffer_t> buffer_list;

    // Target physics posts receives.
    post_receives(buffer_list);

    // Source physics sends buffers to target physics.
    send();

    // Target physics receives buffers from target physics.
    fill_nodes(buffer_list, key);

    // reset the internal communicator
    nemesis::reset_internal_comm();
}


//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Create empty buffers for packing and unpacking data.
 */
void Messenger::calculate_buffer_sizes()
{
    Require (d_map);
    Require ( d_map->partitions().size() == d_map->partition_pts().size() );

    // set the total number of buffers
    d_buffer_sizes.resize(d_map->partitions().size(), 0);

    // Get an iterator into the buffer size vector
    typename Vec_Ord::iterator size_iter = d_buffer_sizes.begin();

    // Get the partition_pts ordinate vector
    const Vec_Ord partition_pts = d_map->partition_pts();

    // Get an iterator into the map's partition pts
    for(typename Vec_Ord::const_iterator iter = partition_pts.begin(),
                                     iter_end = partition_pts.end();
        iter != iter_end; ++iter)
    {
        Check ( size_iter != d_buffer_sizes.end() );

        // Calculate the size and store it
        *size_iter = (*iter) * sizeof(DataType);

        // Increment the size iterator
        ++size_iter;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Post asynchronous receives for the buffers.
 */
void Messenger::post_receives(BufferList &buffer_list)
{
    Check ( d_buffer_sizes.size() == d_map->partitions().size() );

    // Create an iterator over the partitions
    typename Vec_Ord::const_iterator partition_iter = 
        d_map->partitions().begin();

    // Loop over the buffer sizes, create the buffers, and post the receives
    for (typename Vec_Ord::iterator iter = d_buffer_sizes.begin(), 
                                iter_end = d_buffer_sizes.end(); 
         iter != iter_end; ++iter) 
    {
        Check ( partition_iter != d_map->partitions().end() );

        // Buffer source PID
        OrdinateType source = *partition_iter;

        // Create the buffer
        buffer_list.push_back( Message_Buffer_t(source, *iter) );

        // Get the request buffer
        Message_Buffer_t& buffer = buffer_list.back();

        // Post asynchronous receive with this receive buffer
        nemesis::receive_async(buffer.request(), &buffer.buffer()[0], 
                               buffer.buffer().size(), buffer.ordinate());

        // Increment the partition iterators
        ++partition_iter;
    }

    Ensure ( buffer_list.size() == d_map->partitions().size() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do synchronous sends of data with key \a key.
 */
void Messenger::send(const KeyType& key)
{
    typedef typename Vec_Ord::const_iterator        Vec_Ord_Iter;

    // Create an iterator into the buffer sizes
    Vec_Ord_Iter size_iter = d_buffer_sizes.begin();

    // Create an iterator into the map nodes
    typename Vec_Node::const_iterator node_iter = d_map->nodes().begin();
    typename Vec_Node::const_iterator node_iter_end = d_map->nodes().end();

    // make a packer
    denovo::Packer p;

    // Create a buffer for the packer
    Buffer buffer;

    // Loop over the partitions and send
    for (Vec_Ord_Iter partition_iter = d_map->partitions().begin(),
                  partition_iter_end = d_map->partitions().end();
         partition_iter != partition_iter_end; ++partition_iter)
    {
        Check ( size_iter != d_buffer_sizes.end() );

        // resize the buffer
        buffer.resize(*size_iter);

        // buffer destination PID
        OrdinateType destination = *partition_iter;
        Check (destination < nemesis::nodes());

        // register the current buffer with the packer
        p.set_buffer( buffer.size(), &buffer[0] );

        // Pack the data
        while(node_iter->partition() == *partition_iter && 
              node_iter != node_iter_end)
        {
            Check ( node_iter->data_exists(key) );
            p << node_iter->get_data(key);

            ++node_iter;
        }

        // blocking send the buffer to the remote process
        nemesis::send_async( &buffer[0], buffer.size(), destination );

        // Increment the partition and size iterators
        ++size_iter;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Fill the map with the received data.
 */
void Messenger::fill_nodes(BufferList &buffer_list, 
                                       const KeyType& key)
{
    typedef typename BufferList::iterator   BufferList_Iter;
    typedef typename Vec_Node::iterator     Node_Iterator;

    // make a new denovo::Unpacker
    denovo::Unpacker u;

    // Make a piece of data to write into
    DataType data;

    // Keep going until all requests have been processed
    while( !buffer_list.empty() )
    {
        // Find a buffer with a completed communication request
        BufferList_Iter buffer_iter = 
            std::find_if(buffer_list.begin(), buffer_list.end(), 
                         &Message_Buffer_t::complete);

        // If a completed communication request was found, process it
        if( buffer_iter != buffer_list.end() )
        {
            // Get the source partition
            int source = buffer_iter->source();

            // Get the range of nodes that belong to this source partition
            std::pair<Node_Iterator, Node_Iterator> node_range = 
                d_map->find_node_range(source);

            // Get the buffer
            Buffer &buffer = buffer_iter->buffer();

            // Set the buffer for the unpacker
            u.set_buffer( buffer.size(), &buffer[0] );
            
            // Loop over the nodes and insert the information
            for(Node_Iterator node_iter = node_range.first; 
                node_iter != node_range.second; ++node_iter)
            {
                // Unpack the data
                u >> data;

                // Insert into the node
                if( !node_iter->data_exists(key) )
                {
                    node_iter->register_data(key, data);
                }
                else
                {
                    node_iter->get_data(key) = data;
                }
            }

            // Remove this buffer from the list
            buffer_list.erase(buffer_iter);
        }
    }
}

} // end namespace coupler

//---------------------------------------------------------------------------//
//                 end of Messenger.cc
//---------------------------------------------------------------------------//
