//---------------------------------------------------------------------------//
/*! 
 * \file tstCrsMatrixFill.cpp
 * \author Stuart R. Slattery
 * \brief CrsMatrix fill bug tests.
 */
//---------------------------------------------------------------------------//

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CrsMatrix, LocalFill )
{
    typedef double Scalar;
    typedef int LO;
    typedef int GO;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();

    // Make a map.
    int num_local_elements = 1;
    int num_global_elements = comm_size * num_local_elements;
    Teuchos::RCP<const Tpetra::Map<LO,GO> > map =
	Tpetra::createUniformContigMap<LO,GO>( num_global_elements, comm );

    // Make a matrix with a dynamic graph.
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> > A =
	Tpetra::createCrsMatrix<Scalar,LO,GO>( map, num_local_elements );

    // Fill the matrix with global ids.
    Scalar value = 1.0;
    GO global_row = map->getNodeElementList()[0];
    Teuchos::Array<GO> indices(1,global_row);
    Teuchos::Array<Scalar> values(1,value);
    A->insertGlobalValues( global_row, indices(), values() );
    A->fillComplete();

    // Check the matrix.
    Teuchos::ArrayView<const LO> indices_view;
    Teuchos::ArrayView<const Scalar> values_view;
    A->getLocalRowView( 0, indices_view, values_view );
    TEST_EQUALITY( 1, indices_view.size() );
    TEST_EQUALITY( 1, values_view.size() );

    Scalar test_val = value;
    TEST_EQUALITY( test_val, values_view[0] );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CrsMatrix, EveryoneFill )
{
    typedef double Scalar;
    typedef int LO;
    typedef int GO;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();

    // Make a map.
    int num_local_elements = 1;
    int num_global_elements = comm_size * num_local_elements;
    Teuchos::RCP<const Tpetra::Map<LO,GO> > map =
	Tpetra::createUniformContigMap<LO,GO>( num_global_elements, comm );

    // Make a matrix with a dynamic graph.
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> > A =
	Tpetra::createCrsMatrix<Scalar,LO,GO>( map, num_local_elements );

    // Fill the matrix with global ids. Everyone will put something in every
    // row.
    Scalar value = 1.0;
    Teuchos::Array<Scalar> values(1,value);
    Teuchos::Array<GO> indices(1,0);
    for ( GO i = 0; i < num_global_elements; ++i )
    {
	indices[0] = i;
	A->insertGlobalValues( i, indices(), values() );
    }
    A->fillComplete();

    // Check the matrix.
    Teuchos::ArrayView<const LO> indices_view;
    Teuchos::ArrayView<const Scalar> values_view;
    A->getLocalRowView( 0, indices_view, values_view );
    TEST_EQUALITY( 1, indices_view.size() );
    TEST_EQUALITY( 1, values_view.size() );

    Scalar test_val = comm_size * value;
    TEST_EQUALITY( test_val, values_view[0] );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CrsMatrix, NotMyFill )
{
    typedef double Scalar;
    typedef int LO;
    typedef int GO;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
	Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();
    int comm_rank = comm->getRank();

    // Make a map.
    int num_local_elements = 1;
    int num_global_elements = comm_size * num_local_elements;
    Teuchos::RCP<const Tpetra::Map<LO,GO> > map =
	Tpetra::createUniformContigMap<LO,GO>( num_global_elements, comm );

    // Make a matrix with a dynamic graph.
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> > A =
	Tpetra::createCrsMatrix<Scalar,LO,GO>( map, num_local_elements );

    // Fill the matrix with global ids. Fill a row that you do not own.
    Scalar value = 1.0;
    GO global_row = comm_size - comm_rank - 1;
    Teuchos::Array<Scalar> values(1,value);
    Teuchos::Array<GO> indices(1,global_row);
    A->insertGlobalValues( global_row, indices(), values() );
    A->fillComplete();

    // Check the matrix.
    Teuchos::ArrayView<const LO> indices_view;
    Teuchos::ArrayView<const Scalar> values_view;
    A->getLocalRowView( 0, indices_view, values_view );
    TEST_EQUALITY( 1, indices_view.size() );
    TEST_EQUALITY( 1, values_view.size() );

    Scalar test_val = value;
    TEST_EQUALITY( test_val, values_view[0] );
}

//---------------------------------------------------------------------------//
// end tstCrsMatrixFill.cpp
//---------------------------------------------------------------------------//
