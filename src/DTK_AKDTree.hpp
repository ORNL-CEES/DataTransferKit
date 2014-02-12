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
 * \file DTK_AKDTree.hpp
 * \author Stuart R. Slattery
 * \brief AKDTree class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_AKDTREE_HPP
#define DTK_AKDTREE_HPP

#include "DTK_FieldTraits.hpp"
#include "DTK_FieldManager.hpp"
#include "DTK_MeshTraits.hpp"
#include "DTK_MeshManager.hpp"
#include "DTK_BoundingBox.hpp"
#include "DTK_EvaluationPoint.hpp"
#include "DTK_BufferCommunicator.hpp"
#include "DTK_ElementTree.hpp"
#include "DTK_DataBuffer.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class AKDTree 
 * \brief Fully-asynchronous mutually exclusive kd-tree.
 */
//---------------------------------------------------------------------------//
template<class Mesh, class CoordinateField, int DIM>
class AKDTree
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh                                         mesh_type;
    typedef MeshTraits<Mesh>                             MT;
    typedef typename MT::global_ordinal_type             GlobalOrdinal;
    typedef MeshManager<Mesh>                            MeshManagerType;
    typedef Teuchos::RCP<MeshManagerType>                RCP_MeshManager;
    typedef CoordinateField                              coord_field_type;
    typedef FieldManager<CoordinateField>                CoordFieldManagerType;
    typedef Teuchos::RCP<CoordFieldManagerType>          RCP_CoordFieldManager;
    typedef FieldTraits<CoordinateField>                 CFT;
    typedef EvaluationPoint<GlobalOrdinal,DIM>           packet_type;
    typedef typename DataBuffer<packet_type>::BankType   BankType;
    typedef Teuchos::Comm<int>                           Comm;
    typedef Teuchos::CommRequest<int>                    Request;
    //@}

    // Constructor.
    AKDTree( const Teuchos::RCP<const Comm>& comm,
	     const int max_buffer_size,
	     const int buffer_check_frequency,
	     const Teuchos::Array<int>& source_neighbor_ranks,
	     const Teuchos::Array<BoundingBox>& source_neighbor_boxes,
	     const Teuchos::Array<int>& target_neighbor_ranks,
	     const RCP_MeshManager& source_mesh_manager,
	     const RCP_CoordFieldManager& target_coord_manager,
	     const Teuchos::ArrayView<const GlobalOrdinal> target_ordinals );

    // Destructor.
    ~AKDTree() { /* ... */ }

    // Locate the target points in the tree and return the source
    // element/target point pairings in the source decomposition. Target
    // coordinates will be blocked by dimension.
    void locate( Teuchos::Array<GlobalOrdinal>& source_elements,
		 Teuchos::Array<GlobalOrdinal>& target_point_gids,
		 Teuchos::Array<double>& target_coords,
		 const GlobalOrdinal global_num_targets,
		 const double tolerance );

  private:

    // Send a local target point to its source destination.
    void sendTargetPoint( BankType& bank );

    // Process a point in the bank.
    void processBankPoint( BankType& bank,
			   Teuchos::Array<GlobalOrdinal>& source_elements,
			   Teuchos::Array<GlobalOrdinal>& target_point_gids,
			   Teuchos::Array<double>& interleaved_target_coords,
			   const double tolerance );

    // Process incoming messages.
    void processMessages( BankType& bank );

    // Post communications in the binary tree.
    void postTreeCount();

    // Complete outstanding communications in the binary tree at the end of a
    // cycle.
    void completeTreeCount();

    // Update the binary tree count of completed points.
    void updateTreeCount();

    // Send the global finished message to the children.
    void sendCompleteToChildren();

    // Control the termination of a stage.
    void controlTermination();

  private:

    // Master proc enumeration for implementation clarity.
    enum MasterIndicator { MASTER = 0 };

  private:

    // Parallel communicator for this tree.
    Teuchos::RCP<const Comm> d_comm;

    // Parent process.
    int d_parent;

    // Child processes.
    std::pair<int,int> d_children;

    // Buffer communicator.
    BufferCommunicator<packet_type> d_buffer_communicator;

    // Source mesh manager.
    RCP_MeshManager d_source_mesh_manager;

    // Target coordinate manager.
    RCP_CoordFieldManager d_target_coord_manager;

    // Bounding boxes for the source neighbors.
    Teuchos::Array<BoundingBox> d_source_neighbor_boxes;
    
    // Target coordinate global ids.
    Teuchos::ArrayView<const GlobalOrdinal> d_target_ordinals;

    // Master-worker communicator for reporting number of packets.
    Teuchos::RCP<const Comm> d_comm_num_done;

    // Master-worker communicator for reporting completion status.
    Teuchos::RCP<const Comm> d_comm_complete;

    // Master-worker asynchronous communication request handles for number of
    // packets complete.
    std::pair<Teuchos::RCP<Request>,Teuchos::RCP<Request> > d_num_done_handles;

    // Master-worker reports for number of packets complete communications. 
    std::pair<Teuchos::RCP<int>,Teuchos::RCP<int> > d_num_done_report;

    // Request handle for completed work on worker nodes.
    Teuchos::RCP<Request> d_complete_handle;

    // Completion report.
    Teuchos::RCP<int> d_complete_report;

    // Total number of source packets in set.
    GlobalOrdinal d_nh;

    // Total number of packets completed in set.
    Teuchos::RCP<int> d_num_done;
    
    // Number of packets complete in the local problem.
    int d_num_run;

    // Number of target points left to send to their destinations.
    int d_targets_to_send;

    // Boolean-as-integer from completion of search calculation.
    Teuchos::RCP<int> d_complete;

    // Check frequency for data buffer communication.
    int d_check_freq;

    // Element tree for searching the local mesh.
    Teuchos::RCP<ElementTree<Mesh> > d_element_tree;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_AKDTree_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_AKDTREE_HPP

//---------------------------------------------------------------------------//
// end DTK_AKDTree.hpp
//---------------------------------------------------------------------------//

