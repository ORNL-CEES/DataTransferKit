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
 * \file DTK_Distributor.hpp
 * \author Stuart R. Slattery
 * \brief Distributor class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_DISTRIBUTOR_HPP
#define DTK_DISTRIBUTOR_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class Distributor 
 * \brief Data distributor.
 *
 * The Distributor solves the MxN problem and executes the resulting
 * communication plan.
 */
//---------------------------------------------------------------------------//
class Distributor 
{
  public:

    //@{
    //! Typedefs.
    typedef Teuchos::Comm<int>                        Comm;
    typedef Teuchos::CommRequest<int>                 Request;
    //@}

    // Constructor.
    Distributor( const Teuchos::RCP<const Comm>& comm,
		 const int check_frequency = 10 );

    // Destructor.
    ~Distributor() { /* ... */ }

    // Create from sends.
    std::size_t 
    createFromSends( const Teuchos::ArrayView<const int>& export_node_ids );

    // Get the number of messages I will send.
    std::size_t getNumSends() const
    { return d_num_sends; }

    // Get the images to which I am sending.
    Teuchos::ArrayView<const int> getImagesTo() const
    { return d_images_to(); }
    
    // Get the number of packets I am sending.
    Teuchos::ArrayView<const std::size_t> getLengthsTo() const
    { return d_lengths_to(); }

    // Get the number of messages I will receive.
    std::size_t getNumReceives() const
    { return d_num_receives; }

    // Get the images from which I am receiveing.
    Teuchos::ArrayView<const int> getImagesFrom() const
    { return d_images_from(); }
    
    // Get the number of packets I am receiving.
    Teuchos::ArrayView<const std::size_t> getLengthsFrom() const
    { return d_lengths_from(); }

  private:

    // Post communications in the binary tree.
    void postTreeCount();

    // Complete outstanding communications in the binary tree at the end of a
    // cycle.
    void completeTreeCount();

    // Update the binary tree count of completed receives.
    void updateTreeCount();

    // Send the global finished message to the children.
    void sendCompleteToChildren();

    // Control the termination of a stage.
    void controlTermination();

  private:

    // Master proc enumeration for implementation clarity.
    enum MasterIndicator { MASTER = 0 };

  private:

    // Parallel communicator for this set.
    Teuchos::RCP<const Comm> d_comm;

    // Check frequency.
    int d_check_freq;

    // Parent process.
    int d_parent;

    // Child processes.
    std::pair<int,int> d_children;

    // Master-worker communicator tag for reporting number of sends.
    Teuchos::RCP<const Comm> d_comm_num_done;

    // Master-worker communicator tag for reporting completion status.
    Teuchos::RCP<const Comm> d_comm_complete;

    // Master-worker asynchronous communication request handles for number of
    // receives complete.
    std::pair<Teuchos::RCP<Request>,Teuchos::RCP<Request> > d_num_done_handles;

    // Master-worker reports for number of receives complete communications. 
    std::pair<Teuchos::RCP<int>,Teuchos::RCP<int> > d_num_done_report;

    // Request handle for completed work on worker nodes.
    Teuchos::RCP<Request> d_complete_handle;

    // Completion report.
    Teuchos::RCP<int> d_complete_report;

    // Total number of processes that have completed.
    Teuchos::RCP<int> d_num_done;
    
    // Boolean-as-integer from completion of sends.
    Teuchos::RCP<int> d_complete;

    // Total number of messages to be sent.
    int d_global_sends;

    // Number of exports.
    std::size_t d_num_exports;

    // Number of processes I am sending to.
    std::size_t d_num_sends;

    // List of process Ids to which I am sending including myself.
    Teuchos::Array<int> d_images_to;

    // Starting index of the block of packets to send to each process.
    Teuchos::Array<std::size_t> d_starts_to;

    // Length (in number of packets) I will send to each process.
    Teuchos::Array<std::size_t> d_lengths_to;

    // The number of messages I will receive.
    std::size_t d_num_receives;

    // List of process Ids from which I am receiving including myself.
    Teuchos::Array<int> d_images_from;

    // Starting index of the block of packets received from each process.
    Teuchos::Array<std::size_t> d_starts_from;

    // Length (in number of packets) I will receive from each process.
    Teuchos::Array<std::size_t> d_lengths_from;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_DISTRIBUTOR_HPP

//---------------------------------------------------------------------------//
// end DTK_Distributor.hpp
//---------------------------------------------------------------------------//

