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
 * \file   tstMovingLeastSquare.cpp
 * \author Stuart R. Slattery
 * \brief  MovingLeastSquare tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <DTK_MeshFreeInterpolator.hpp>
#include <DTK_MeshFreeInterpolatorFactory.hpp>

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
// Floating point epsilon.
//---------------------------------------------------------------------------//

const double epsilon = 0.4;

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( MovingLeastSquare, dim_1_test )
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    int rank = comm->getRank();
    int size = comm->getSize();
    int inverse_rank = size - rank - 1;

    int dim = 1;
    int num_src_points = 2;
    int num_src_coords = dim*num_src_points;

    // Source coordinates.
    Teuchos::Array<double> src_coords(num_src_coords);

    // Point 1.
    src_coords[0] = 2.0*rank;

    // Point 2.
    src_coords[1] = 2.0*rank + 1.0;

    // Source Data.
    Teuchos::Array<double> src_data( num_src_points );
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_data[i] = 3.0*rank + 1.0;
    }

    int num_tgt_points = 2;
    int num_tgt_coords = dim*num_tgt_points;
    Teuchos::Array<double> tgt_coords( num_tgt_coords );
    tgt_coords[0] = 2.0*inverse_rank + 0.4;
    tgt_coords[1] = 2.0*inverse_rank + 0.6;

    double radius = 1.0;

    Teuchos::RCP<DataTransferKit::MeshFreeInterpolator> interpolator = 
	DataTransferKit::MeshFreeInterpolatorFactory::create<int>(
	    comm, "Moving Least Square", "Wendland", 2, dim );

    interpolator->setProblem( src_coords(), tgt_coords(), radius );

    Teuchos::Array<double> tgt_data( num_tgt_points );
    interpolator->interpolate( src_data(), 1, src_data.size(),
			       tgt_data(), 1, tgt_data.size() );

    for ( int i = 0; i < num_tgt_points; ++i )
    {
	TEST_FLOATING_EQUALITY( tgt_data[i], 3.0*rank + 1.0, epsilon );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( MovingLeastSquare, dim_2_test )
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    int rank = comm->getRank();
    int size = comm->getSize();
    int inverse_rank = size - rank - 1;

    int dim = 2;
    int num_src_points = 4;
    int num_src_coords = dim*num_src_points;

    // Source coordinates.
    Teuchos::Array<double> src_coords(num_src_coords);

    // Point 1.
    src_coords[0] = 2.0*rank;
    src_coords[1] = 0.0;

    // Point 2.
    src_coords[2] = 2.0*rank + 1.0;
    src_coords[3] = 0.0;

    // Point 3.
    src_coords[4] = 2.0*rank + 1.0;
    src_coords[5] = 1.0;

    // Point 4.
    src_coords[6] = 2.0*rank;
    src_coords[7] = 1.0;

    // Source Data.
    Teuchos::Array<double> src_data( num_src_points );
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_data[i] = 3.0*rank + 1.0;
    }

    int num_tgt_points = 2;
    int num_tgt_coords = dim*num_tgt_points;
    Teuchos::Array<double> tgt_coords( num_tgt_coords );
    tgt_coords[0] = 2.0*inverse_rank + 0.4;
    tgt_coords[1] = 0.5;
    tgt_coords[2] = 2.0*inverse_rank + 0.6;
    tgt_coords[3] = 0.5;

    double radius = 1.0;

    Teuchos::RCP<DataTransferKit::MeshFreeInterpolator> interpolator = 
	DataTransferKit::MeshFreeInterpolatorFactory::create<int>(
	    comm, "Moving Least Square", "Wendland", 2, dim );

    interpolator->setProblem( src_coords(), tgt_coords(), radius );

    Teuchos::Array<double> tgt_data( num_tgt_points );
    interpolator->interpolate( src_data(), 1, src_data.size(),
			       tgt_data(), 1, tgt_data.size() );
    for ( int i = 0; i < num_tgt_points; ++i )
    {
	TEST_FLOATING_EQUALITY( tgt_data[i], 3.0*rank + 1.0, epsilon );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( MovingLeastSquare, dim_3_test )
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    int rank = comm->getRank();
    int size = comm->getSize();
    int inverse_rank = size - rank - 1;

    int dim = 3;
    int num_src_points = 8;
    int num_src_coords = dim*num_src_points;

    // Source coordinates.
    Teuchos::Array<double> src_coords(num_src_coords);

    // Point 1.
    src_coords[0] = 2.0*rank;
    src_coords[1] = 0.0;
    src_coords[2] = 0.0;

    // Point 2.
    src_coords[3] = 2.0*rank + 1.0;
    src_coords[4] = 0.0;
    src_coords[5] = 0.0;

    // Point 3.
    src_coords[6] = 2.0*rank + 1.0;
    src_coords[7] = 1.0;
    src_coords[8] = 0.0;

    // Point 4.
    src_coords[9] = 2.0*rank;
    src_coords[10] = 1.0;
    src_coords[11] = 0.0;

    // Point 5.
    src_coords[12] = 2.0*rank;
    src_coords[13] = 0.0;
    src_coords[14] = 1.0;

    // Point 6.
    src_coords[15] = 2.0*rank + 1.0;
    src_coords[16] = 0.0;
    src_coords[17] = 1.0;

    // Point 7.
    src_coords[18] = 2.0*rank + 1.0;
    src_coords[19] = 1.0;
    src_coords[20] = 1.0;

    // Point 8.
    src_coords[21] = 2.0*rank;
    src_coords[22] = 1.0;
    src_coords[23] = 1.0;

    // Source Data.
    Teuchos::Array<double> src_data( num_src_points );
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_data[i] = 3.0*rank + 1.0;
    }

    int num_tgt_points = 2;
    int num_tgt_coords = dim*num_tgt_points;
    Teuchos::Array<double> tgt_coords( num_tgt_coords );
    tgt_coords[0] = 2.0*inverse_rank + 0.4;
    tgt_coords[1] = 0.5;
    tgt_coords[2] = 0.5;
    tgt_coords[3] = 2.0*inverse_rank + 0.6;
    tgt_coords[4] = 0.5;
    tgt_coords[5] = 0.5;

    double radius = 1.0;

    Teuchos::RCP<DataTransferKit::MeshFreeInterpolator> interpolator = 
	DataTransferKit::MeshFreeInterpolatorFactory::create<int>(
	    comm, "Moving Least Square", "Wendland", 2, dim );

    interpolator->setProblem( src_coords(), tgt_coords(), radius );

    Teuchos::Array<double> tgt_data( num_tgt_points );
    interpolator->interpolate( src_data(), 1, src_data.size(),
			       tgt_data(), 1, tgt_data.size() );
    for ( int i = 0; i < num_tgt_points; ++i )
    {
	TEST_FLOATING_EQUALITY( tgt_data[i], 3.0*rank + 1.0, epsilon );
    }
}

//---------------------------------------------------------------------------//
// end tstMovingLeastSquare.cpp
//---------------------------------------------------------------------------//
