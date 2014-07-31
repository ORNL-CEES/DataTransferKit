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
 * \file   tstOrthogonalPolynomialMLSProblem.cpp
 * \author Stuart R. Slattery
 * \brief  OrthogonalPolynomialMLSProblem tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <DTK_MovingLeastSquare.hpp>
#include <DTK_WuBasis.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_ParameterList.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// Get the default communicator.
Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<int>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<int>() );
#endif
}

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( OrthogonalPolynomialMLSProblem, 3by3_test )
{
    // Setup a 3x3 grid.
    const int dim = 3;
    int num_source = 9;
    Teuchos::Array<double> source_centers( num_source*dim );
    Teuchos::Array<unsigned> source_ids( num_source );
    for ( int i = 0; i < num_source; ++i )
    {
	source_ids[i] = i;
    }
    
    // P1.
    source_centers[0] = -1.0;
    source_centers[1] = -1.0;
    source_centers[2] = 0.0;

    // P2.
    source_centers[3] = 0.0;
    source_centers[4] = -1.0;
    source_centers[5] = 0.0;

    // P3.
    source_centers[6] = 1.0;
    source_centers[7] = -1.0;
    source_centers[8] = 0.0;

    // P4.
    source_centers[9] = -1.0;
    source_centers[10] = 0.0;
    source_centers[11] = 0.0;

    // P5.
    source_centers[12] = 0.0;
    source_centers[13] = 0.0;
    source_centers[14] = 0.0;

    // P6.
    source_centers[15] = 1.0;
    source_centers[16] = 0.0;
    source_centers[17] = 0.0;

    // P7.
    source_centers[18] = -1.0;
    source_centers[19] = 1.0;
    source_centers[20] = 0.0;

    // P8.
    source_centers[21] = 0.0;
    source_centers[22] = 1.0;
    source_centers[23] = 0.0;

    // P9.
    source_centers[24] = 1.0;
    source_centers[25] = 1.0;
    source_centers[26] = 0.0;

    // Evaluate the basis at the center.
    Teuchos::Array<double> target_center( dim, 0.0 );

    // Build a moving least square interpolation.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    DataTransferKit::MovingLeastSquare<DataTransferKit::WuBasis<4>,int,dim>
	mls( comm, "OP" );
    mls.setProblem( source_centers(), target_center(), 1.5 );

    Teuchos::Array<double> source_data( num_source, 1.0 );
    Teuchos::Array<double> target_data( 1, 0.0 );
    mls.interpolate( source_data(), target_data(), 1 );
    TEST_EQUALITY( target_data[0], 1.0 );
}

//---------------------------------------------------------------------------//
// end tstOrthogonalPolynomialMLSProblem.cpp
//---------------------------------------------------------------------------//
