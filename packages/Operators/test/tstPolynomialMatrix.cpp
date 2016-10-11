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
 * \file   tstPolynomialMatrix.cpp
 * \author Stuart R. Slattery
 * \brief  Basis function unit tests.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <DTK_PolynomialMatrix.hpp>

#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// Get the default communicator.
template <class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal>> getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp( new Teuchos::SerialComm<Ordinal>() );
#endif
}

//---------------------------------------------------------------------------//
// Floating point epsilon.
//---------------------------------------------------------------------------//

const double epsilon_abs = 1.0e-12;
const double epsilon_rel = 1.0e-7 / 100.0; // percent tolerance divided by 100

//---------------------------------------------------------------------------//
// Unit tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( PolynomialMatrix, polynomial_matrix_apply )
{
    // Make an equivalent polynomial matrix and CrsMatrix and apply it to a
    // multivector.
    Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm<int>();
    int local_size = 100;
    int global_size = local_size * comm->getSize();
    int num_vec = 3;
    int poly_size = 10;

    // Create a random polynomial.
    Teuchos::RCP<const Tpetra::Map<int, DataTransferKit::SupportId>> row_map =
        Tpetra::createUniformContigMap<int, DataTransferKit::SupportId>(
            global_size, comm );
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        P = Tpetra::createMultiVector<double, int, DataTransferKit::SupportId>(
            row_map, poly_size );
    P->randomize();

    // Create the CrsMatrix version of the polynomial.
    Teuchos::RCP<const Tpetra::Map<int, DataTransferKit::SupportId>> col_map =
        Tpetra::createLocalMap<int, DataTransferKit::SupportId>( poly_size,
                                                                 comm );
    Teuchos::RCP<Tpetra::CrsMatrix<double, int, DataTransferKit::SupportId>>
        P_crs = Teuchos::rcp(
            new Tpetra::CrsMatrix<double, int, DataTransferKit::SupportId>(
                row_map, col_map, poly_size, Tpetra::StaticProfile ) );
    Teuchos::Array<int> crs_indices( poly_size );
    for ( int j = 0; j < poly_size; ++j )
    {
        crs_indices[j] = j;
    }
    Teuchos::Array<double> crs_values( poly_size );
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double>> poly_view =
        P->get2dView();
    for ( int i = 0; i < local_size; ++i )
    {
        for ( int j = 0; j < poly_size; ++j )
        {
            crs_values[j] = poly_view[j][i];
        }
        P_crs->insertLocalValues( i, crs_indices(), crs_values() );
    }
    P_crs->fillComplete();

    // Create the PolynomialMatrix version of the polynomial.
    DataTransferKit::PolynomialMatrix P_poly_mat( P, row_map, row_map );

    // Build a random vector to apply the matrices to.
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        X = Tpetra::createMultiVector<double, int, DataTransferKit::SupportId>(
            row_map, num_vec );
    X->randomize();

    // Apply the CrsMatrix.
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        Y_crs =
            Tpetra::createMultiVector<double, int, DataTransferKit::SupportId>(
                row_map, num_vec );
    Y_crs->randomize();
    P_crs->apply( *X, *Y_crs, Teuchos::NO_TRANS );

    // Apply the polynomial matrix.
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        Y_poly_mat =
            Tpetra::createMultiVector<double, int, DataTransferKit::SupportId>(
                row_map, num_vec );
    Y_poly_mat->randomize();
    P_poly_mat.apply( *X, *Y_poly_mat, Teuchos::NO_TRANS );

    // Compare the results.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double>> y_crs_view =
        Y_crs->get2dView();
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double>> y_pm_view =
        Y_poly_mat->get2dView();
    for ( int i = 0; i < num_vec; ++i )
    {
        for ( int j = 0; j < local_size; ++j )
        {
            TEST_FLOATING_EQUALITY( y_crs_view[i][j], y_pm_view[i][j],
                                    epsilon_rel );
            TEST_COMPARE( std::abs( y_crs_view[i][j] - y_pm_view[i][j] ), <=,
                          epsilon_abs );
        }
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( PolynomialMatrix, polynomial_matrix_transpose_apply )
{
    // Make an equivalent polynomial matrix and CrsMatrix and apply it to a
    // multivector.
    Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm<int>();
    int local_size = 100;
    int global_size = local_size * comm->getSize();
    int num_vec = 3;
    int poly_size = 10;

    // Create a random polynomial.
    Teuchos::RCP<const Tpetra::Map<int, DataTransferKit::SupportId>> row_map =
        Tpetra::createUniformContigMap<int, DataTransferKit::SupportId>(
            global_size, comm );
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        P = Tpetra::createMultiVector<double, int, DataTransferKit::SupportId>(
            row_map, poly_size );
    P->randomize();

    // Create the CrsMatrix version of the polynomial.
    Teuchos::RCP<const Tpetra::Map<int, DataTransferKit::SupportId>> col_map =
        Tpetra::createLocalMap<int, DataTransferKit::SupportId>( poly_size,
                                                                 comm );
    Teuchos::RCP<Tpetra::CrsMatrix<double, int, DataTransferKit::SupportId>>
        P_crs = Teuchos::rcp(
            new Tpetra::CrsMatrix<double, int, DataTransferKit::SupportId>(
                row_map, col_map, poly_size, Tpetra::StaticProfile ) );
    Teuchos::Array<int> crs_indices( poly_size );
    for ( int j = 0; j < poly_size; ++j )
    {
        crs_indices[j] = j;
    }
    Teuchos::Array<double> crs_values( poly_size );
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double>> poly_view =
        P->get2dView();
    for ( int i = 0; i < local_size; ++i )
    {
        for ( int j = 0; j < poly_size; ++j )
        {
            crs_values[j] = poly_view[j][i];
        }
        P_crs->insertLocalValues( i, crs_indices(), crs_values() );
    }
    P_crs->fillComplete();

    // Create the PolynomialMatrix version of the polynomial.
    DataTransferKit::PolynomialMatrix P_poly_mat( P, row_map, row_map );

    // Build a random vector to apply the matrices to.
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        X = Tpetra::createMultiVector<double, int, DataTransferKit::SupportId>(
            row_map, num_vec );
    X->randomize();

    // Transpose apply the CrsMatrix.
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        Y_crs =
            Tpetra::createMultiVector<double, int, DataTransferKit::SupportId>(
                row_map, num_vec );
    Y_crs->randomize();
    P_crs->apply( *X, *Y_crs, Teuchos::TRANS );

    // Transpose apply the polynomial matrix.
    Teuchos::RCP<Tpetra::MultiVector<double, int, DataTransferKit::SupportId>>
        Y_poly_mat =
            Tpetra::createMultiVector<double, int, DataTransferKit::SupportId>(
                row_map, num_vec );
    Y_poly_mat->randomize();
    P_poly_mat.apply( *X, *Y_poly_mat, Teuchos::TRANS );

    // Check the results.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double>> y_crs_view =
        Y_crs->get2dView();
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double>> y_pm_view =
        Y_poly_mat->get2dView();
    for ( int i = 0; i < num_vec; ++i )
    {
        for ( int j = 0; j < local_size; ++j )
        {
            TEST_FLOATING_EQUALITY( y_crs_view[i][j], y_pm_view[i][j],
                                    epsilon_rel );
            TEST_COMPARE( std::abs( y_crs_view[i][j] - y_pm_view[i][j] ), <=,
                          epsilon_abs );
        }
    }
}

//---------------------------------------------------------------------------//
// end tstPolynomialMatrix.cpp
//---------------------------------------------------------------------------//
