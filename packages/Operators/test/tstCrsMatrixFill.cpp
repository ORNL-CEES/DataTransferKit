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

#include <DTK_Types.hpp>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CrsMatrix, LocalFill )
{
    typedef double Scalar;
    typedef int LO;
    typedef DataTransferKit::SupportId GO;

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
    typedef DataTransferKit::SupportId GO;

    // Get the communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
        Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();

    // Make a map.
    int num_local_elements = 1;
    GO num_global_elements = comm_size * num_local_elements;
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
    typedef DataTransferKit::SupportId GO;

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
