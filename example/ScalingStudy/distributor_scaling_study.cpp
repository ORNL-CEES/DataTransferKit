//---------------------------------------------------------------------------//
/*!
 * \file distributor_scaling_study.cpp
 * \author Stuart R. Slattery
 * \brief Distributor scaling study.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <set>

#include <DTK_Distributor.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Tpetra_Distributor.hpp"

//---------------------------------------------------------------------------//
// Weak scaling study driver.
//---------------------------------------------------------------------------//
int main(int argc, char* argv[])
{
    // Setup communication.
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int my_rank = comm->getRank();
    int my_size = comm->getSize();

    // Decide who we are sending to.
    int num_sends = std::atoi( argv[1] );
    Teuchos::Array<int> send_ranks( num_sends );
    srand( my_rank + my_rank*237 );
    for ( int n = 0; n < num_sends; ++n )
    {
	rand();
	send_ranks[n] = std::floor( my_size * double(rand())/double(RAND_MAX) );
    }

    // Create the tpetra distributor.
    Tpetra::Distributor tpetra_distributor( comm );
    int num_import_tpetra = 0;
    {
	Teuchos::RCP<Teuchos::Time> tpetra_timer = 
	    Teuchos::TimeMonitor::getNewCounter("Tpetra");
	Teuchos::TimeMonitor monitor(*tpetra_timer);
	num_import_tpetra = tpetra_distributor.createFromSends( send_ranks() );
    }

    // Create the DTK distributor.
    DataTransferKit::Distributor dtk_distributor( comm, 8 );
    int num_import_dtk = 0;
    {
	Teuchos::RCP<Teuchos::Time> dtk_timer = 
	    Teuchos::TimeMonitor::getNewCounter("DataTransferKit");
	Teuchos::TimeMonitor monitor(*dtk_timer);
	num_import_dtk = dtk_distributor.createFromSends( send_ranks() );
    }

    // Check DTK against Tpetra.
    int local_test_fail = 0;

    // Check the number of sends and receives.
    int tpetra_self_message = tpetra_distributor.hasSelfMessage() ? 1 : 0;
    if ( num_import_tpetra != num_import_dtk ) ++local_test_fail;
    if ( tpetra_distributor.getNumSends() + tpetra_self_message !=
		   dtk_distributor.getNumSends() ) ++local_test_fail;
    if ( tpetra_distributor.getNumReceives() + tpetra_self_message !=
		   dtk_distributor.getNumReceives() ) ++local_test_fail;

    // Check the images we are sending to.
    Teuchos::Array<int> tpetra_images_to( tpetra_distributor.getImagesTo() );
    std::sort( tpetra_images_to.begin(), tpetra_images_to.end() );
    Teuchos::Array<int> dtk_images_to( dtk_distributor.getImagesTo() );
    std::sort( dtk_images_to.begin(), dtk_images_to.end() );
    int num_to_images = dtk_images_to.size();
    for ( int i = 0; i < num_to_images; ++i )
    {
	if ( tpetra_images_to[i] != dtk_images_to[i] ) ++local_test_fail;
    }

    // Check the images we are receiving from.
    Teuchos::Array<int> tpetra_images_from( tpetra_distributor.getImagesFrom() );
    std::sort( tpetra_images_from.begin(), tpetra_images_from.end() );
    Teuchos::Array<int> dtk_images_from( dtk_distributor.getImagesFrom() );
    std::sort( dtk_images_from.begin(), dtk_images_from.end() );
    int num_from_images = dtk_images_from.size();
    for ( int i = 0; i < num_from_images; ++i )
    {
	if ( tpetra_images_from[i] != dtk_images_from[i] ) ++local_test_fail;
    }

    // Check the number of packets we are sending to each image.
    Teuchos::Array<std::size_t> tpetra_lengths_to( tpetra_distributor.getLengthsTo() );
    std::sort( tpetra_lengths_to.begin(), tpetra_lengths_to.end() );
    Teuchos::Array<std::size_t> dtk_lengths_to( dtk_distributor.getLengthsTo() );
    std::sort( dtk_lengths_to.begin(), dtk_lengths_to.end() );
    int num_to_lengths = dtk_lengths_to.size();
    for ( int i = 0; i < num_to_lengths; ++i )
    {
	if ( tpetra_lengths_to[i] != dtk_lengths_to[i] ) ++local_test_fail;
    }

    // Check the number of packets we are receiving from each image.
    Teuchos::Array<std::size_t> tpetra_lengths_from( tpetra_distributor.getLengthsFrom() );
    std::sort( tpetra_lengths_from.begin(), tpetra_lengths_from.end() );
    Teuchos::Array<std::size_t> dtk_lengths_from( dtk_distributor.getLengthsFrom() );
    std::sort( dtk_lengths_from.begin(), dtk_lengths_from.end() );
    int num_from_lengths = dtk_lengths_from.size();
    for ( int i = 0; i < num_from_lengths; ++i )
    {
	if ( tpetra_lengths_from[i] != dtk_lengths_from[i] ) ++local_test_fail;
    }

    // Compute the global number of test failures.
    int global_test_fail = 0;
    Teuchos::reduceAll( *comm,  Teuchos::REDUCE_SUM,
			local_test_fail, Teuchos::Ptr<int>(&global_test_fail) );

    // Output results
    if ( my_rank == 0 )
    {
    	std::cout << "==================================================" 
    		  << std::endl;
    	std::cout << "DTK Distributor weak scaling study" << std::endl;
    	std::cout << "Number of processors: " << my_size << std::endl;
    	std::cout << "Number of sends:      " << num_sends << std::endl;
	std::cout << "Test failures:        " << global_test_fail << std::endl;
    	std::cout << "==================================================" 
    		  << std::endl;
    }
    comm->barrier();

    Teuchos::TableFormat &format = Teuchos::TimeMonitor::format();
    format.setPrecision(5);
    Teuchos::TimeMonitor::summarize();

    return 0;
}

//---------------------------------------------------------------------------//
// end distributor_scaling_study.cpp
//---------------------------------------------------------------------------//

