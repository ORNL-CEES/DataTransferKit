//---------------------------------------------------------------------------//
/*!
 * \file   tstMpiTagConsistency.cpp
 * \author Stuart Slattery
 * \brief  Mpi Tags unit tests.
 */
//---------------------------------------------------------------------------//

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
// This is a 2 processor test.
TEUCHOS_UNIT_TEST( CommTools, mpi_tag_consistency )
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm_default = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm_default->getSize();
    int comm_rank = comm_default->getRank();

    // Split the main communicator into 2 overlapping groups
    Teuchos::Array<int> sub_ranks_1(2);
    sub_ranks_1[0] = 0;
    sub_ranks_1[1] = 1;
    Teuchos::RCP<const Teuchos::Comm<int> > comm_1 = 
	comm_default->createSubcommunicator( sub_ranks_1() );

    Teuchos::Array<int> sub_ranks_2(1);
    sub_ranks_2[0] = 0;
    Teuchos::RCP<const Teuchos::Comm<int> > comm_2 = 
	comm_default->createSubcommunicator( sub_ranks_2() );

    // Create another communicator.
    Teuchos::Array<int> sub_ranks_3(comm_size);
    for ( int i = 0; i < comm_size; ++i )
    {
	sub_ranks_3[i] = i;
    }
    Teuchos::RCP<const Teuchos::Comm<int> > comm_3 =
	comm_default->createSubcommunicator( sub_ranks_3() );

    // Get my mpi tag for comm 3.
    int my_tag = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >( 
        comm_3 )->getTag();

    // Collect the tags for comm 3.
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

    // Test the tags.
    TEST_EQUALITY( tag1, tag2 );
}

//---------------------------------------------------------------------------//
// end tstMpiTagConsistency.cpp
//---------------------------------------------------------------------------//
