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
 * \file   parallel_search.cpp
 * \author Stuart Slattery
 * \brief  STK parallel search example.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <cstdlib>

#include "DTK_STKMeshManager.hpp"
#include "DTK_STKMeshEntityPredicates.hpp"
#include "DTK_ParallelSearch.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <Tpetra_MultiVector.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_topology/topology.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

//---------------------------------------------------------------------------//
// Example driver.
//---------------------------------------------------------------------------//
int main(int argc, char* argv[])
{
    // INITIALIZATION
    // --------------

    // Setup communication.
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();

    // Read in command line options.
    std::string xml_input_filename;
    Teuchos::CommandLineProcessor clp(false);
    clp.setOption( "xml-in-file",
		   &xml_input_filename,
		   "The XML file to read into a parameter list" );
    clp.parse(argc,argv);

    // Build the parameter list from the xml input.
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile(
	xml_input_filename, Teuchos::inoutArg(*plist) );

    // Read command-line options
    std::string source_mesh_input_file = 
	plist->get<std::string>("Source Mesh Input File");
    std::string source_mesh_part_name = 
	plist->get<std::string>("Source Mesh Part");
    std::string target_mesh_input_file = 
	plist->get<std::string>("Target Mesh Input File");
    std::string target_mesh_part_name = 
	plist->get<std::string>("Target Mesh Part");

    // Get the raw mpi communicator (basic typedef in STK).
    Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm = 
	Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( comm );
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm = 
	mpi_comm->getRawMpiComm();
    stk::ParallelMachine parallel_machine = (*opaque_comm)();

    
    // SOURCE MESH READ
    // ----------------

    // Load the source mesh.
    stk::io::StkMeshIoBroker src_broker( parallel_machine );
    std::size_t src_input_index = src_broker.add_mesh_database(
	source_mesh_input_file, "exodus", stk::io::READ_MESH );
    src_broker.set_active_mesh( src_input_index );
    src_broker.create_input_mesh();

    // Create the source bulk data.
    src_broker.populate_bulk_data();
    Teuchos::RCP<stk::mesh::BulkData> src_bulk_data = 
	Teuchos::rcpFromRef( src_broker.bulk_data() );

    // Create a selector from the source part.
    stk::mesh::Part* src_part = 
	src_broker.meta_data().get_part( source_mesh_part_name );
    stk::mesh::Selector src_stk_selector( *src_part );

    
    // TARGET MESH READ
    // ----------------

    // Load the target mesh.
    stk::io::StkMeshIoBroker tgt_broker( parallel_machine );
    std::size_t tgt_input_index = tgt_broker.add_mesh_database(
    	target_mesh_input_file, "exodus", stk::io::READ_MESH );
    tgt_broker.set_active_mesh( tgt_input_index );
    tgt_broker.create_input_mesh();

    // Create the target bulk data.
    tgt_broker.populate_bulk_data();
    Teuchos::RCP<stk::mesh::BulkData> tgt_bulk_data = 
    	Teuchos::rcpFromRef( tgt_broker.bulk_data() );

    // Create a selector from the target part.
    stk::mesh::Part* tgt_part = 
	tgt_broker.meta_data().get_part( target_mesh_part_name );
    stk::mesh::Selector tgt_stk_selector( *tgt_part );    

    
    // PARALLEL SEARCH
    // -----------------------

    // Get the search parameters.
    Teuchos::ParameterList& dtk_params = plist->sublist( "DataTransferKit" );
    Teuchos::ParameterList& search_params = dtk_params.sublist( "Search" );
    
    // Create a manager for the source mesh.
    Teuchos::RCP<DataTransferKit::ClientManager> src_manager =
	Teuchos::rcp( new DataTransferKit::STKMeshManager(src_bulk_data) );

    // Create a manager for the target mesh.
    Teuchos::RCP<DataTransferKit::ClientManager> tgt_manager =
	Teuchos::rcp( new DataTransferKit::STKMeshManager(tgt_bulk_data) );

    // Create a parallel search over the elements in the source mesh.
    DataTransferKit::PredicateFunction src_select_predicate =
	DataTransferKit::STKSelectorPredicate( src_stk_selector );
    DataTransferKit::EntityIterator src_element_iterator =
	src_manager->entitySet()->entityIterator( 3, src_select_predicate );
    DataTransferKit::ParallelSearch parallel_search( comm,
						     3,
						     src_element_iterator,
						     src_manager->localMap(),
						     search_params );

    // Search for the location of the target mesh nodes in the source mesh
    // elements.
    DataTransferKit::PredicateFunction tgt_select_predicate =
	DataTransferKit::STKSelectorPredicate( tgt_stk_selector );
    DataTransferKit::EntityIterator tgt_node_iterator =
	tgt_manager->entitySet()->entityIterator( 0, tgt_select_predicate );
    parallel_search.search( tgt_node_iterator,
			    tgt_manager->localMap(),
			    search_params );


    // PARALLEL SEARCH OUTPUT
    // ----------------------
    
    // In the parallel decomposition of the source mesh, print out which
    // target nodes were found in which elements and the parallel ranks that
    // own those target points.
    Teuchos::Array<DataTransferKit::EntityId> tgt_ids;
    Teuchos::ArrayView<const double> tgt_param_coords;
    DataTransferKit::EntityIterator src_elems_begin =
	src_element_iterator.begin();
    DataTransferKit::EntityIterator src_elems_end =
	src_element_iterator.end();
    for ( auto src_elem = src_elems_begin;
	  src_elem != src_elems_end;
	  ++src_elem )
    {
	std::cout << "SOURCE MESH ELEMENT " << src_elem->id() << std::endl;
	
	// Get the ids of the target mesh nodes found in this source mesh
	// element.
	parallel_search.getRangeEntitiesFromDomain( src_elem->id(), tgt_ids );

	// Get the target node search data in the source element.
	int num_tgt = tgt_ids.size();
	for ( int i = 0; i < num_tgt; ++i )
	{
	    // Get the parametric coordinates of the target node in the source
	    // element. 
	    parallel_search.rangeParametricCoordinatesInDomain(
		src_elem->id(), tgt_ids[i], tgt_param_coords );
	    
	    std::cout << "TARGET NODE " << tgt_ids[i] << std::endl;
	    std::cout << "TARGET RANK "
		      << parallel_search.rangeEntityOwnerRank( tgt_ids[i] )
		      << std::endl;
	    std::cout << "TARGET PARAMETRIC COORDS IN ELEMENT "
		      << tgt_param_coords;
	}

	std::cout << std::endl;
    }

    // In the parallel decomposition of the target mesh, print out the source
    // elements in which the target nodes were found and the owning parallel
    // ranks of the source elements.
    Teuchos::Array<DataTransferKit::EntityId> src_ids;
    Teuchos::Array<int> src_ranks;
    DataTransferKit::EntityIterator tgt_nodes_begin =
	tgt_node_iterator.begin();
    DataTransferKit::EntityIterator tgt_nodes_end =
	tgt_node_iterator.end();
    for ( auto tgt_node = tgt_nodes_begin;
	  tgt_node != tgt_nodes_end;
	  ++tgt_node )
    {
	// Get the ids of the source mesh elements in which this target node
	// was found.
	parallel_search.getDomainEntitiesFromRange( tgt_node->id(), src_ids );

	// Get the parallel ranks that own the source elements.
	int num_src = src_ids.size();
	src_ranks.resize( num_src );
	for ( int i = 0; i < num_src; ++i )
	{
	    src_ranks[i] = parallel_search.domainEntityOwnerRank( src_ids[i] );
	}

	// Print the target node ids and ranks.
	std::cout << "TARGET NODE " << tgt_node->id() << std::endl;
	std::cout << "SOURCE ELEMENT IDS " << src_ids << std::endl;
	std::cout << "SOURCE RANKS       " << src_ranks << std::endl;
	std::cout << std::endl;	
    }

    // Print the ids of the target mesh nodes on this process that were not
    // found in any source mesh element.
    std::cout << "MISSED TARGET NODES "
	      << parallel_search.getMissedRangeEntityIds()
	      << std::endl << std::endl;
}

//---------------------------------------------------------------------------//
// end tstSTK_Mesh.cpp
//---------------------------------------------------------------------------//
