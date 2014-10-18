//---------------------------------------------------------------------------//
/*!
 * \file tstDistributor.cpp
 * \author Stuart R. Slattery
 * \brief Tpetra Distributor unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Tpetra_Distributor.hpp>

//---------------------------------------------------------------------------//
// MPI Setup
//---------------------------------------------------------------------------//

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
// Tests
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( Distributor, distributor_test )
{
    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    Teuchos::Array<int> export_data( my_size );
    Teuchos::Array<int> export_procs( my_size );
    for ( int i = 0; i < my_size; ++i )
    {
	export_data[i] = i;
	export_procs[i] = i;
    }

    Tpetra::Distributor distributor( getDefaultComm<int>() );

    int num_import = distributor.createFromSends( export_procs() );
    TEST_ASSERT( num_import == my_size );

    Teuchos::ArrayView<const int> export_data_view = export_data();

    Teuchos::Array<int> import_data( num_import );
    Teuchos::ArrayView<int> import_data_view = import_data();

    distributor.doPostsAndWaits( export_data_view, 1, import_data_view );

    for ( int i = 0; i < num_import; ++i )
    {
	TEST_ASSERT( import_data[i] == my_rank );
    }
}

//---------------------------------------------------------------------------//
// end tstDistributor.cpp
//---------------------------------------------------------------------------//

