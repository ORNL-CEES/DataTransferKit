//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tstDistributor.cpp
 * \author Stuart Slattery
 * \date   Wed May 25 12:36:14 2011
 * \brief  Distributor class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <cstdlib>
#include <algorithm>

#include "DTK_Distributor.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "Tpetra_Distributor.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// Get the default communicator.
template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( Distributor, distributor_test )
{
    // Comm setup.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm<int>();
    int size = comm->getSize();
    int rank = comm->getRank();

    // Make a series of random send ranks.
    int num_sends = 1000;
    Teuchos::Array<int> export_ids( num_sends );
    srand( rank*num_sends );
    for ( int i = 0; i < num_sends; ++i )
    {
	export_ids[i] = std::floor(
	    size * double(rand())/double(RAND_MAX) );
    }

    // Construct the communication plan with Tpetra.
    Tpetra::Distributor tpetra_distributor( comm );
    int num_import_tpetra = tpetra_distributor.createFromSends( export_ids() );

    // Construct the communication plan with DTK.
    DataTransferKit::Distributor dtk_distributor( comm );
    int num_import_dtk = dtk_distributor.createFromSends( export_ids() );
    

    // Check DTK against Tpetra.
    int tpetra_self_message = tpetra_distributor.hasSelfMessage() ? 1 : 0;
    TEST_EQUALITY( num_import_tpetra, num_import_dtk );
    TEST_EQUALITY( tpetra_distributor.getNumSends() + tpetra_self_message,
		   dtk_distributor.getNumSends() );
    TEST_EQUALITY( tpetra_distributor.getNumReceives() + tpetra_self_message,
		   dtk_distributor.getNumReceives() );

    Teuchos::Array<int> tpetra_images_to( tpetra_distributor.getImagesTo() );
    std::sort( tpetra_images_to.begin(), tpetra_images_to.end() );
    Teuchos::Array<int> dtk_images_to( dtk_distributor.getImagesTo() );
    std::sort( dtk_images_to.begin(), dtk_images_to.end() );
    int num_to_images = dtk_images_to.size();
    for ( int i = 0; i < num_to_images; ++i )
    {
	TEST_EQUALITY( tpetra_images_to[i], dtk_images_to[i] );
    }

    Teuchos::Array<int> tpetra_images_from( tpetra_distributor.getImagesFrom() );
    std::sort( tpetra_images_from.begin(), tpetra_images_from.end() );
    Teuchos::Array<int> dtk_images_from( dtk_distributor.getImagesFrom() );
    std::sort( dtk_images_from.begin(), dtk_images_from.end() );
    int num_from_images = dtk_images_from.size();
    for ( int i = 0; i < num_from_images; ++i )
    {
	TEST_EQUALITY( tpetra_images_from[i], dtk_images_from[i] );
    }

    Teuchos::Array<std::size_t> tpetra_lengths_to( tpetra_distributor.getLengthsTo() );
    std::sort( tpetra_lengths_to.begin(), tpetra_lengths_to.end() );
    Teuchos::Array<std::size_t> dtk_lengths_to( dtk_distributor.getLengthsTo() );
    std::sort( dtk_lengths_to.begin(), dtk_lengths_to.end() );
    int num_to_lengths = dtk_lengths_to.size();
    for ( int i = 0; i < num_to_lengths; ++i )
    {
	TEST_EQUALITY( tpetra_lengths_to[i], dtk_lengths_to[i] );
    }

    Teuchos::Array<std::size_t> tpetra_lengths_from( tpetra_distributor.getLengthsFrom() );
    std::sort( tpetra_lengths_from.begin(), tpetra_lengths_from.end() );
    Teuchos::Array<std::size_t> dtk_lengths_from( dtk_distributor.getLengthsFrom() );
    std::sort( dtk_lengths_from.begin(), dtk_lengths_from.end() );
    int num_from_lengths = dtk_lengths_from.size();
    for ( int i = 0; i < num_from_lengths; ++i )
    {
	TEST_EQUALITY( tpetra_lengths_from[i], dtk_lengths_from[i] );
    }
}

//---------------------------------------------------------------------------//
//                        end of tstDistributor.cpp
//---------------------------------------------------------------------------//
