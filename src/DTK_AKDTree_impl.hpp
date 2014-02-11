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
 * \file DTK_AKDTree_impl.hpp
 * \author Stuart R. Slattery
 * \brief AKDTree class implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_AKDTREE_IMPL_HPP
#define DTK_AKDTREE_IMPL_HPP

#include "DTK_DBC.hpp"
#include "DTK_CommTools.hpp"
#include "DTK_FieldTools.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_OrdinalTraits.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Mesh, class CoordinateField, int DIM>
AKDTree<Mesh,CoordinateField,DIM>::AKDTree( 
    const Teuchos::RCP<const Comm>& comm,
    const int max_buffer_size,
    const int buffer_check_frequency,
    const Teuchos::Array<int>& source_neighbor_ranks,
    const Teuchos::Array<BoundingBox>& source_neighbor_boxes,
    const Teuchos::Array<int>& target_neighbor_ranks,
    const RCP_MeshManager& source_mesh_manager,
    const RCP_CoordFieldManager& target_coord_manager,
    const Teuchos::ArrayView<const GlobalOrdinal> target_ordinals )
    : d_comm( comm )
    , d_parent( Teuchos::OrdinalTraits<int>::invalid() )
    , d_children( Teuchos::OrdinalTraits<int>::invalid(),
                  Teuchos::OrdinalTraits<int>::invalid() )
    , d_buffer_communicator( source_neighbor_ranks, target_neighbor_ranks,
			     d_comm, max_buffer_size )
    , d_source_mesh_manager( source_mesh_manager )
    , d_target_coord_manager( target_coord_manager )
    , d_source_neighbor_boxes( source_neighbor_boxes )
    , d_target_ordinals( target_ordinals )
    , d_num_done_report( Teuchos::rcp(new int(0)), Teuchos::rcp(new int(0)) )
    , d_complete_report( Teuchos::rcp(new int(0)) )
    , d_num_done( Teuchos::rcp(new int(0)) )
    , d_complete( Teuchos::rcp(new int(0)) )
    , d_check_freq( buffer_check_frequency )
{
    DTK_REQUIRE( Teuchos::nonnull(d_comm) );

    // Set the duplicate communicators. This is how we get around not having
    // access to message tags through the abstract Teuchos::Comm interface. We
    // are constructing a separate messaging space for each of these
    // bookeeping operations.
    d_comm_num_done = d_comm->duplicate();
    d_comm_complete = d_comm->duplicate();

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

    // Build an element tree for searching the local mesh.
    d_source_mesh_manager->buildIndexing();
    d_element_tree = Teuchos::rcp( 
	new ElementTree<Mesh>(d_source_mesh_manager) );

    DTK_ENSURE( Teuchos::nonnull(d_comm_num_done) );
    DTK_ENSURE( Teuchos::nonnull(d_comm_complete) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Locate the target points in the tree and return the source
 * element/target point pairings in the source decomposition. Target
 * coordinates will be blocked by dimension.
 */
template<class Mesh, class CoordinateField, int DIM>
void AKDTree<Mesh,CoordinateField,DIM>::locate( 
    Teuchos::Array<GlobalOrdinal>& source_elements,
    Teuchos::Array<GlobalOrdinal>& target_point_gids,
    Teuchos::Array<double>& target_coords,
    const double tolerance )
{
    // Initialize.
    *d_complete = 0;
    *d_num_done = 0;
    d_num_run = 0;
    d_targets_to_send = d_target_ordinals.size();

    // Get the number of local target points.
    Teuchos::reduceAll<int,GlobalOrdinal>( 
	*d_comm,
	Teuchos::REDUCE_SUM,
	d_targets_to_send,
	Teuchos::Ptr<GlobalOrdinal>(&d_nh) );

    // Make a temporary array of interleaved target coordinates.
    Teuchos::Array<double> interleaved_target_coords;

    // Barrier before continuing.
    d_comm->barrier();

    // Create a data bank.
    BankType bank;
    DTK_CHECK( bank.empty() );

    // Everyone posts receives for data buffers to get started.
    d_buffer_communicator.post();

    // Post asynchronous communcations in the binary tree for packet counts.
    postTreeCount();

    // Process points in the tree until completion.
    while ( !(*d_complete) )
    {
	// Process local target points.
	if ( 0 < d_targets_to_send )
	{
	    sendTargetPoint( bank );
	    --d_targets_to_send;
	    ++d_num_run;
	}

	// If the source is empty, process the local target points in the bank.
	else if ( !bank.empty() )
	{
	    processBankPoint( bank, 
			      source_elements, 
			      target_point_gids, 
			      interleaved_target_coords, 
			      tolerance );
	    ++d_num_run;
	}

	// If we're out of source and bank target points or have hit the check
	// frequency, process incoming messages.
	if ( ((0 == d_targets_to_send) && bank.empty()) ||  
	     d_num_run == d_check_freq )
	{
	    processMessages( bank );
	    d_num_run = 0;
	}

	// If everything looks like it is finished locally, report through
	// the tree to check if transport is done.
	if ( (0 == d_targets_to_send) && bank.empty() )
	{
	    controlTermination();
	}
    }

    // Barrier before continuing.
    d_comm->barrier();

    // Complete the binary tree outstanding communication.
    completeTreeCount();

    // End all communication and free all buffers.
    DTK_CHECK( !d_buffer_communicator.sendBufferSize() );
    DTK_CHECK( bank.empty() );
    d_buffer_communicator.end();
    DTK_CHECK( !d_buffer_communicator.sendStatus() );
    DTK_CHECK( !d_buffer_communicator.receiveStatus() );

    // Reorder the target coordinates to be block by dimension.
    GlobalOrdinal num_targets = target_point_gids.size();
    target_coords.resize( interleaved_target_coords.size() );
    for ( GlobalOrdinal i = 0; i < num_targets; ++i )
    {
	for ( int d = 0; d < DIM; ++d )
	{
	    target_coords[ num_targets*d + i ] =
		interleaved_target_coords[ DIM*i + d ];
	}
    }

    // Barrier before completion.
    d_comm->barrier();
    DTK_ENSURE( 0 == d_targets_to_send );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send a target point to its source destination.
 */
template<class Mesh, class CoordinateField, int DIM>
void AKDTree<Mesh,CoordinateField,DIM>::sendTargetPoint( BankType& bank )
{
    DTK_REQUIRE( 0 < d_targets_to_send );
    DTK_REQUIRE( Teuchos::nonnull(d_target_coord_manager) );

    // Get a point from the source.
    GlobalOrdinal num_targets = d_target_ordinals.size();
    GlobalOrdinal point_id = num_targets - d_targets_to_send;
    Teuchos::Array<double> coords(DIM);
    Teuchos::ArrayRCP<const double> coords_view = 
	FieldTools<CoordinateField>::view( *d_target_coord_manager->field() );
    for ( int d = 0; d < DIM; ++d )
    {
	coords[d] = coords_view[ num_targets*d + point_id ];
    }
    Teuchos::RCP<packet_type> point = Teuchos::rcp(
	new packet_type(d_target_ordinals[point_id], coords()) );

    // Send it to the neighboring boxes in which it resides.
    bool point_in_domain = false;
    for ( unsigned n = 0; n < d_source_neighbor_boxes.size(); ++n )
    {
	if ( d_source_neighbor_boxes[n].pointInBox(coords) )
	{
	    point_in_domain = true;
	    d_buffer_communicator.communicate( point, n );
	}
    }

    // If the point was outside the domain, update the completion count.
    if ( !point_in_domain )
    {
	++(*d_num_done);
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process a point in the bank.
 */
template<class Mesh, class CoordinateField, int DIM>
void AKDTree<Mesh,CoordinateField,DIM>::processBankPoint(
    BankType& bank,
    Teuchos::Array<GlobalOrdinal>& source_elements,
    Teuchos::Array<GlobalOrdinal>& target_point_gids,
    Teuchos::Array<double>& interleaved_target_coords,
    const double tolerance )
{
    DTK_REQUIRE( !bank.empty() );
    DTK_REQUIRE( Teuchos::nonnull(d_source_mesh_manager) );
    DTK_REQUIRE( Teuchos::nonnull(d_element_tree) );

    // Get a point from the bank.
    Teuchos::RCP<packet_type> point = bank.top();
    bank.pop();
    DTK_CHECK( Teuchos::nonnull(point) );

    // Search in the element tree for the point.
    GlobalOrdinal element_gid = 
	Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
    if ( d_element_tree->findPoint(point->coords(),element_gid,tolerance) )
    {
	source_elements.push_back( element_gid );
	target_point_gids.push_back( point->gid() );
	for ( int d = 0; d < DIM; ++d )
	{
	    interleaved_target_coords.push_back( point->coordinate(d) );
	}
    }

    // Update the completion count.
    ++(*d_num_done);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Process incoming messages.
 */
template<class Mesh, class CoordinateField, int DIM>
void AKDTree<Mesh,CoordinateField,DIM>::processMessages( BankType& bank )
{
    // Check for incoming packets.
    d_buffer_communicator.checkAndPost(bank);

    // Add to the packet completed tally 
    updateTreeCount();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Post communications in the binary tree.
 */
template<class Mesh, class CoordinateField, int DIM>
void AKDTree<Mesh,CoordinateField,DIM>::postTreeCount()
{
    // Post a receive from the first child for packet count data.
    if ( d_children.first != Teuchos::OrdinalTraits<int>::invalid() )
    {
	*d_num_done_report.first = 0;
	d_num_done_handles.first = Teuchos::ireceive<int,int>(
	    *d_comm_num_done, d_num_done_report.first, d_children.first );
    }

    // Post a receive from the second child for packet count data.
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
template<class Mesh, class CoordinateField, int DIM>
void AKDTree<Mesh,CoordinateField,DIM>::completeTreeCount()
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
 * \brief Update the binary tree count of completed packets.
 */
template<class Mesh, class CoordinateField, int DIM>
void AKDTree<Mesh,CoordinateField,DIM>::updateTreeCount()
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
            DTK_CHECK( *d_num_done <= d_nh );

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
            DTK_CHECK( *d_num_done <= d_nh );

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
template<class Mesh, class CoordinateField, int DIM>
void AKDTree<Mesh,CoordinateField,DIM>::sendCompleteToChildren()
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
template<class Mesh, class CoordinateField, int DIM>
void AKDTree<Mesh,CoordinateField,DIM>::controlTermination()
{
    // Send any partially full buffers.
    d_buffer_communicator.send();

    // Update the packet count from the children.
    updateTreeCount();

    // MASTER checks for completion.
    if ( d_comm->getRank() == MASTER ) 
    {
        if ( *d_num_done == d_nh )
        {
            *d_complete = 1;
            sendCompleteToChildren();
        }
    }

    // Other nodes send the number of packets completed to parent and check
    // to see if a message has arrived from the parent indicating that
    // transport has been completed.
    else
    {
        // Send completed number of packets to parent.
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

#endif // end DTK_AKDTREE_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_AKDTree_impl.hpp
//---------------------------------------------------------------------------//

