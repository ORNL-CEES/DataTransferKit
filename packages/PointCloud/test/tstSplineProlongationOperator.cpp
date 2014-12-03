//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \file   tstSplineProlongationOperator.cpp
 * \author Stuart R. Slattery
 * \brief  Basis function unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <limits>

#include <DTK_SplineProlongationOperator.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

//---------------------------------------------------------------------------//
// Unit tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( SplineProlongationOperator, polynomial_matrix_apply )
{
    // Make an equivalent polynomial matrix and CrsMatrix and apply it to a
    // multivector.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int local_size = 100;
    int global_size = local_size * comm->getSize();
    int num_vec = 3;
    int offset = (comm->getRank() == 0) ? 4 : 0;

    // Create a map.
    Teuchos::RCP<const Tpetra::Map<int,int> > map = 
	Tpetra::createUniformContigMap<int,int>( global_size, comm );

    // Create a prolongator.
    DataTransferKit::SplineProlongationOperator<double,int> 
	prolongation_op( offset, map );

    // Check the prolongator.
    TEST_EQUALITY( prolongation_op.getRangeMap()->getNodeNumElements(),
		   Teuchos::as<unsigned>(offset + local_size) );

    // Build a random vector to prolongate.
    Teuchos::RCP<Tpetra::MultiVector<double,int,int> > X =
	Tpetra::createMultiVector<double,int,int>( map, num_vec );
    X->randomize();

    // Build a prolongated vector.
    Teuchos::RCP<Tpetra::MultiVector<double,int,int> > Y =
	Tpetra::createMultiVector<double,int,int>( 
	    prolongation_op.getRangeMap(), num_vec );
    Y->putScalar( 1.0 );

    // Prolongate.
    prolongation_op.apply( *X, *Y );

    // Compare the results.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double> > X_view = X->get2dView();
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double> > Y_view = Y->get2dView();
    for ( int i = 0; i < num_vec; ++i )
    {
	for ( int j = 0; j < offset; ++j )
	{
	    TEST_EQUALITY( Y_view[i][j], 0.0 );
	}
	for ( int j = offset; j < local_size; ++j )
	{
	    TEST_EQUALITY( Y_view[i][j], X_view[i][j-offset] );
	}
    }
}

//---------------------------------------------------------------------------//
// end tstSplineProlongationOperator.cpp
//---------------------------------------------------------------------------//
