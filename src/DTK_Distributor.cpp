//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \file DTK_Distributor.cpp
 * \author Stuart R. Slattery
 * \brief Distributor class implementation.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "DTK_Distributor.hpp"
#include "DTK_DBC.hpp"
#include "DTK_CommTools.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#include <mpi.h>
#endif

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Distributor::Distributor( const Teuchos::RCP<const Comm>& comm )
    : d_comm( comm )
    , d_parent( Teuchos::OrdinalTraits<int>::invalid() )
    , d_children( Teuchos::OrdinalTraits<int>::invalid(),
                  Teuchos::OrdinalTraits<int>::invalid() )
    , d_num_done_report( Teuchos::rcp(new int(0)), Teuchos::rcp(new int(0)) )
    , d_complete_report( Teuchos::rcp(new int(0)) )
    , d_num_done( Teuchos::rcp(new int(0)) )
    , d_complete( Teuchos::rcp(new int(0)) )
{
    DTK_REQUIRE( !d_comm.is_null() );

    // We are constructing a separate messaging space for each of these
    // bookeeping operations. Right now we are using the copy constructor to
    // simply change the tag of the messages. If isend and ireceive had single
    // object semantics with a tag option in Teuchos::CommHelpers we would use
    // that instead.
#ifdef HAVE_MPI
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm = 
	Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( d_comm );
    Teuchos::RCP< const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm = 
	mpi_comm->getRawMpiComm();
    d_comm_num_done = Teuchos::rcp( new Teuchos::MpiComm<int>(opaque_comm) );
    d_comm_complete = Teuchos::rcp( new Teuchos::MpiComm<int>(opaque_comm) );
#else
    d_comm_num_done = comm;
    d_comm_complete = comm;
#endif

    // Get the comm parameters.
    int my_rank = d_comm->getRank();
    int my_size = d_comm->getSize();

    // Get the parent process. MASTER has no parent.
    if ( my_rank != MASTER )
    {
	if ( my_rank % 2 == 0 )
	{
	    d_parent = ( my_rank / 2 ) - 1;
	}
	else
	{
	    d_parent = ( my_rank - 1 ) / 2;
	}
    }
	 
    // Get the first child process.
    int child_1 = ( my_rank * 2 ) + 1;
    if ( child_1 < my_size )
    {
	d_children.first = child_1;
    }

    // Get the second child process.
    int child_2 = child_1 + 1;
    if ( child_2 < my_size )
    {
	d_children.second = child_2;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Create the communication plan from sends.
 */
std::size_t Distributor::createFromSends(
    const Teuchos::ArrayView<const int>& export_node_ids )
{
    // Determine which ranks we are sending to and how many packets we are
    // sending.
    d_num_exports = export_node_ids.size();
    d_images_to = Teuchos::Array<int>( export_node_ids );
    std::sort( d_images_to.begin(), d_images_to.end() );
    Teuchos::Array<int>::iterator image_it = 
	std::unique( d_images_to.begin(), d_images_to.end() );
    d_images_to.resize( std::distance(d_images_to.begin(),image_it) );
    d_num_sends = d_images_to.size();
    d_lengths_to.resize( d_num_sends );
    Teuchos::Array<std::size_t>::iterator length_it;
    for ( image_it = d_images_to.begin(),
	 length_it = d_lengths_to.begin();
	  image_it != d_images_to.end(); 
	  ++image_it, ++length_it )
    {
	*length_it = std::count( export_node_ids.begin(),
				 export_node_ids.end(),
				 *image_it );
    }

    // Reduce the total number of sends to the root rank so we can determine
    // the termination condition.
#ifdef HAVE_MPI
    Teuchos::RCP< const Teuchos::MpiComm<int> > mpi_comm = 
	Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm<int> >( d_comm );
    Teuchos::RCP< const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm = 
	mpi_comm->getRawMpiComm();
    MPI_Comm raw_comm = (*opaque_comm)();
    d_global_sends = 0;
    MPI_Reduce( 
	&d_num_sends, &d_global_sends, 1, MPI_INT, MPI_SUM, MASTER, raw_comm );
#else
    d_global_sends = d_num_sends;
#endif

    // Everyone posts sends to get started.
    int sends_not_done = 0;
    bool local_sends_done = false;
    Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest<int> > > 
	send_requests( d_num_sends );
    for ( std::size_t m = 0; m < d_num_sends; ++m )
    {
	send_requests[m] = Teuchos::isend<int,std::size_t>( 
	    *d_comm, 
	    Teuchos::RCP<std::size_t>(&d_lengths_to[m],false), 
	    d_images_to[m] );
    }

    // Post a receive.
    int any_rank = 0;
#ifdef HAVE_MPI
    any_rank = MPI_ANY_SOURCE;
#endif
    int send_rank = -1;
    Teuchos::RCP<std::size_t> receive_packet = 
	Teuchos::rcp( new std::size_t(0) );
    Teuchos::RCP<Teuchos::CommRequest<int> > receive_request =
	Teuchos::ireceive<int,std::size_t>( *d_comm, receive_packet, any_rank );

    // Post asynchronous communications in the binary tree for completion count.
    postTreeCount();

    // Keep checking for receives until all send requests have been completed.
    *d_complete = 0;
    *d_num_done = 0;
    while ( !(*d_complete) )
    {
	// Check the receive status. If we got something, process it,
	// repost, and update the completion count.
	if ( CommTools::isRequestCompleteWithRank(receive_request,send_rank) )
	{
	    // Get the rank the message came from.
	    d_images_from.push_back( send_rank );
	    send_rank = -1;
	    DTK_CHECK( d_images_from.back() >= 0 );
	    DTK_CHECK( d_images_from.back() < d_comm->getSize() );

	    // Get the number of packets.
	    d_lengths_from.push_back(*receive_packet);
	    DTK_CHECK( d_lengths_from.back() > 0 );

	    // Repost the request.
	    receive_request = Teuchos::ireceive<int,std::size_t>( 
		*d_comm, receive_packet, any_rank );

	    // Update the completion count.
	    *d_num_done += 1;
	}

	// Check to see if the local sends have been completed.
	if ( !local_sends_done )
	{
	    sends_not_done = 0;
	    for ( std::size_t m = 0; m < d_num_sends; ++m )
	    {
		if ( !CommTools::isRequestComplete(send_requests[m]) )
		{
		    ++sends_not_done;
		}
	    }

	    if ( !local_sends_done && 0 == sends_not_done )
	    {
		local_sends_done = true;
	    }
	}

	// If all of our sends have been received then check for the
	// termination condition.
	if ( local_sends_done )
	{
	    controlTermination();
	}
    }

    // Cancel all oustanding send requests.
    Teuchos::waitAll( *d_comm, send_requests() );

    // Cancel the receive request.
    receive_request = Teuchos::null;

    // End the binary tree outstanding communication.
    completeTreeCount();

    // Get the number of receives.
    d_num_receives = d_images_from.size();

    // Return the number of imports this process will receive.
    std::size_t num_imports = 
	std::accumulate( d_lengths_from.begin(), d_lengths_from.end(), 0 );
    return num_imports;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Post communications in the binary tree.
 */
void Distributor::postTreeCount()
{
    // Post a receive from the first child for send count data.
    if ( d_children.first != Teuchos::OrdinalTraits<int>::invalid() )
    {
	*d_num_done_report.first = 0;
	d_num_done_handles.first = Teuchos::ireceive<int,int>(
	    *d_comm_num_done, d_num_done_report.first, d_children.first );
    }

    // Post a receive from the second child for send count data.
    if ( d_children.second != Teuchos::OrdinalTraits<int>::invalid() )
    {
	*d_num_done_report.second = 0;
	d_num_done_handles.second = Teuchos::ireceive<int,int>(
	    *d_comm_num_done, d_num_done_report.second, d_children.second );
    }

    // Post a receive from parent for transport completion.
    if ( d_parent != Teuchos::OrdinalTraits<int>::invalid() )
    {
	d_complete_handle = Teuchos::ireceive<int,int>( 
	    *d_comm_complete, d_complete_report, d_parent );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Complete outstanding communications in the binary tree at the end of
 * a cycle.
 */
void Distributor::completeTreeCount()
{
    Teuchos::Ptr<Teuchos::RCP<Request> > request_ptr;

    // Children nodes send the finish message to the parent.
    if ( d_parent != Teuchos::OrdinalTraits<int>::invalid() )
    {
	Teuchos::RCP<int> clear = Teuchos::rcp( new int(1) );
	Teuchos::RCP<Request> finish = Teuchos::isend<int,int>(
	    *d_comm_num_done, clear, d_parent );
	request_ptr = 
	    Teuchos::Ptr<Teuchos::RCP<Request> >(&finish);
	Teuchos::wait( *d_comm_num_done, request_ptr );
	DTK_CHECK( finish.is_null() );
    }

    // Parent will wait for first child node to clear communication.
    if ( d_children.first != Teuchos::OrdinalTraits<int>::invalid() )
    {
        request_ptr = 
            Teuchos::Ptr<Teuchos::RCP<Request> >(&d_num_done_handles.first);
        Teuchos::wait( *d_comm_num_done, request_ptr );
        DTK_CHECK( d_num_done_handles.first.is_null() );
    }

    // Parent will wait for second child node to clear communication.
    if ( d_children.second != Teuchos::OrdinalTraits<int>::invalid() )
    {
        request_ptr = 
            Teuchos::Ptr<Teuchos::RCP<Request> >(&d_num_done_handles.second);
        Teuchos::wait( *d_comm_num_done, request_ptr );
        DTK_CHECK( d_num_done_handles.second.is_null() );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Update the binary tree count of completed receives.
 */
void Distributor::updateTreeCount()
{
    // Check for received reports of updated counts from first child.
    if ( d_children.first != Teuchos::OrdinalTraits<int>::invalid() )
    {
        // Receive completed reports and repost.
        if ( CommTools::isRequestComplete(d_num_done_handles.first) )
        {
            DTK_CHECK( *(d_num_done_report.first) > 0 );
            d_num_done_handles.first = Teuchos::null;

            // Add to the running total.
            *d_num_done += *(d_num_done_report.first);

            // Repost.
            d_num_done_handles.first = Teuchos::ireceive<int,int>(
                *d_comm_num_done, d_num_done_report.first, d_children.first );
        }
    }

    // Check for received reports of updated counts from second child.
    if ( d_children.second != Teuchos::OrdinalTraits<int>::invalid() )
    {
        // Receive completed reports and repost.
        if ( CommTools::isRequestComplete(d_num_done_handles.second) )
        {
            DTK_CHECK( *(d_num_done_report.second) > 0 );
            d_num_done_handles.second = Teuchos::null;

            // Add to the running total.
            *d_num_done += *(d_num_done_report.second);

            // Repost.
            d_num_done_handles.second = Teuchos::ireceive<int,int>(
                *d_comm_num_done, d_num_done_report.second, d_children.second );
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send the global finished message to the children.
 */
void Distributor::sendCompleteToChildren()
{
    // Child 1
    if ( d_children.first != Teuchos::OrdinalTraits<int>::invalid() )
    {
        Teuchos::RCP<Request> complete = Teuchos::isend<int,int>(
            *d_comm_complete, d_complete, d_children.first );
        Teuchos::Ptr<Teuchos::RCP<Request> > request_ptr(&complete);
        Teuchos::wait( *d_comm_complete, request_ptr );
        DTK_CHECK( complete.is_null() );
    }

    // Child 2
    if ( d_children.second != Teuchos::OrdinalTraits<int>::invalid() )
    {
        Teuchos::RCP<Request> complete = Teuchos::isend<int,int>(
            *d_comm_complete, d_complete, d_children.second );
        Teuchos::Ptr<Teuchos::RCP<Request> > request_ptr(&complete);
        Teuchos::wait( *d_comm_complete, request_ptr );
        DTK_CHECK( complete.is_null() );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Control the termination of a stage.
 */
void Distributor::controlTermination()
{
    // Update the send count from the children.
    updateTreeCount();

    // MASTER checks for completion. We are done if the number of receives
    // equals the number of sends.
    if ( d_comm->getRank() == MASTER ) 
    {
	DTK_CHECK( *d_num_done <= d_global_sends );
        if ( *d_num_done == d_global_sends )
        {
            *d_complete = 1;
            sendCompleteToChildren();
        }
    }

    // Other nodes send the number of ranks that have completed to parent and
    // check to see if a message has arrived from the parent indicating that
    // the process has been completed.
    else
    {
        // Send completed number of receives to parent.
        if ( *d_num_done > 0 )
        {
            Teuchos::RCP<Request> report = Teuchos::isend<int,int>(
                *d_comm_num_done, d_num_done, d_parent );
            Teuchos::Ptr<Teuchos::RCP<Request> > request_ptr(&report);
            Teuchos::wait( *d_comm_num_done, request_ptr );
            DTK_CHECK( report.is_null() );
            *d_num_done = 0;
        } 

        // Check for completion status from parent.
        if ( CommTools::isRequestComplete(d_complete_handle) )
        {
            DTK_CHECK( *d_complete_report ==  1 );
            d_complete_handle = Teuchos::null;
            *d_complete = 1;
            sendCompleteToChildren();
        }
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_Distributor.cpp
//---------------------------------------------------------------------------//

