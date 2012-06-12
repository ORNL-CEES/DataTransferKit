//---------------------------------------------------------------------------//
/*!
 * \file tstTpetraExport.cpp
 * \author Stuart R. Slattery
 * \brief Tpetra TpetraExport unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_MultiVector.hpp>

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

TEUCHOS_UNIT_TEST( Export, tpetra_export_test )
{
    /* 
       Export Map.
       Rank 0:      [ 0 1 2 3 ]
       Rank 1:            [ 3 4 5 6 ]
       Rank 2:                  [ 6 7 8 9 ]
       Rank 3:                        [ 9 10 11 12 ]
       
       Import Map.
       Rank 0:                        [ 9 10 11 12 ]
       Rank 1:                  [ 6 7 8 9 ]
       Rank 2:            [ 3 4 5 6 ]
       Rank 3:      [ 0 1 2 3 ]
     */

    int my_rank = getDefaultComm<int>()->getRank();
    int my_size = getDefaultComm<int>()->getSize();

    if ( my_size == 4 )
    {
	Teuchos::Array<int> export_ids(4), import_ids(4);
	if ( my_rank == 0 )
	{
	    export_ids[0] = 0; import_ids[0] = 9;
	    export_ids[1] = 1; import_ids[1] = 10;
	    export_ids[2] = 2; import_ids[2] = 11;
	    export_ids[3] = 3; import_ids[3] = 12;
	}
	if ( my_rank == 1 )
	{
	    export_ids[0] = 3; import_ids[0] = 6;
	    export_ids[1] = 4; import_ids[1] = 7;
	    export_ids[2] = 5; import_ids[2] = 8;
	    export_ids[3] = 6; import_ids[3] = 9;
	}
	if ( my_rank == 2 )
	{
	    export_ids[0] = 6; import_ids[0] = 3;
	    export_ids[1] = 7; import_ids[1] = 4;
	    export_ids[2] = 8; import_ids[2] = 5;
	    export_ids[3] = 9; import_ids[3] = 6;
	}
	if ( my_rank == 3 )
	{
	    export_ids[0] = 9;  import_ids[0] = 0;
	    export_ids[1] = 10; import_ids[1] = 1;
	    export_ids[2] = 11; import_ids[2] = 2;
	    export_ids[3] = 12; import_ids[3] = 3;
	}
	getDefaultComm<int>()->barrier();

	Teuchos::ArrayView<const int> export_view = export_ids();
	Teuchos::RCP< const Tpetra::Map<int> > export_map =
	    Tpetra::createNonContigMap<int>( 
		export_view, getDefaultComm<int>() );

	Teuchos::ArrayView<const int> import_view = import_ids();
	Teuchos::RCP< const Tpetra::Map<int> > import_map =
	    Tpetra::createNonContigMap<int>( 
		import_view, getDefaultComm<int>() );

	Tpetra::Import<int> importer( export_map, import_map );

	Tpetra::MultiVector<int> export_vector( export_map, 1 );
	export_vector.putScalar( 2 );

	Tpetra::MultiVector<int> import_vector( import_map, 1 );
	import_vector.doImport( export_vector, importer, Tpetra::INSERT );

	Teuchos::ArrayRCP<const int> export_vector_view = 
	    export_vector.get1dView();
	Teuchos::ArrayRCP<const int>::const_iterator export_iterator;
	for ( int i = 0; i < my_size; ++i )
	{
	    if ( my_rank == i )
	    {
		std::cout << "EXPORT RANK " << my_rank << ": ";
		for ( export_iterator = export_vector_view.begin();
		      export_iterator != export_vector_view.end();
		      ++export_iterator )
		{
		    std::cout << *export_iterator << " ";
		}
		std::cout << std::endl;
	    }
	    getDefaultComm<int>()->barrier();
	}

	Teuchos::ArrayRCP<const int> import_vector_view = 
	    import_vector.get1dView();
	Teuchos::ArrayRCP<const int>::const_iterator import_iterator;
	for ( int i = 0; i < my_size; ++i )
	{
	    if ( my_rank == i )
	    {
		std::cout << "IMPORT RANK " << my_rank << ": ";
		for ( import_iterator = import_vector_view.begin();
		      import_iterator != import_vector_view.end();
		      ++import_iterator )
		{
		    std::cout << *import_iterator << " ";
		}
		std::cout << std::endl;
	    }
	    getDefaultComm<int>()->barrier();
	}
    }
}

//---------------------------------------------------------------------------//
// end tstTpetraExport.cpp
//---------------------------------------------------------------------------//

