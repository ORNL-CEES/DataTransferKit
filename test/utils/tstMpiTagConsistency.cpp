//---------------------------------------------------------------------------//
/*!
 * \file   tstCommTools.cpp
 * \author Stuart Slattery
 * \brief  CommTools class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <DTK_CommTools.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_Ptr.hpp"

#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

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
// Tests.
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( CommTools, mpi_tag_consistency_1 )
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm_default = 
	Teuchos::DefaultComm<int>::getComm();

    // Split the main communicator into 2 separate groups
    Teuchos::Array<int> sub_ranks_1, sub_ranks_2;
    int comm_size = comm_default->getSize();
    int comm_rank = comm_default->getRank();

    TEST_EQUALITY( comm_size, 2 );

    for ( int n = 0; n < comm_size; ++n )
    {
        if ( n % 2 == 0 )
        {
            sub_ranks_1.push_back(n);
        }

        sub_ranks_2.push_back(n);
    }

    // Generate the 1 and 2 communicators from the sub ranks.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_1 = 
	comm_default->createSubcommunicator( sub_ranks_1() );
    Teuchos::RCP<const Teuchos::Comm<int> > comm_2 = 
	comm_default->createSubcommunicator( sub_ranks_2() );

    // Unite the communicators.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_union;
    DataTransferKit::CommTools::unite( comm_1, comm_2, comm_union );

    // Get my mpi tag.
    int my_tag = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( 
        comm_union )->getTag();

    // Collect the tags.
    int tag1 = 0;
    if ( comm_rank == 0 )
    {
        tag1 = my_tag;
    }
    comm_default->barrier();
    Teuchos::broadcast( *comm_default, 0, Teuchos::Ptr<int>(&tag1) );

    int tag2 = 0;
    if ( comm_rank == 1 )
    {
        tag2 = my_tag;
    }
    comm_default->barrier();
    Teuchos::broadcast( *comm_default, 1, Teuchos::Ptr<int>(&tag2) );

    TEST_EQUALITY( tag1, tag2 );
}

//---------------------------------------------------------------------------//
// end tstCommTools.cpp
//---------------------------------------------------------------------------//
