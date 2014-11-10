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
 * \file   tstCenterDistributor.cpp
 * \author Stuart R. Slattery
 * \brief  Center distributor tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <DTK_CenterDistributor.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

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
TEUCHOS_UNIT_TEST( CenterDistributor, dim_2_test )
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    int rank = comm->getRank();
    int size = comm->getSize();
    int inverse_rank = size - rank - 1;

    int dim = 2;
    int num_src_points = 10;
    int num_src_coords = dim*num_src_points;

    Teuchos::Array<double> src_coords(num_src_coords);
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_coords[dim*i] = 1.0*i;
	src_coords[dim*i+1] = 2.0*rank;
    }

    int num_tgt_points = 2;
    int num_tgt_coords = dim*num_tgt_points;
    Teuchos::Array<double> tgt_coords( num_tgt_coords );
    tgt_coords[0] = 4.9;
    tgt_coords[1] = 2.0*inverse_rank;
    tgt_coords[2] = 11.4;
    tgt_coords[3] = 2.0*inverse_rank;

    double radius = 1.5;

    Teuchos::Array<double> tgt_decomp_src;

    DataTransferKit::CenterDistributor<2> distributor( 
	comm, src_coords(), tgt_coords(), radius, tgt_decomp_src );

    int num_import = 6;
    TEST_EQUALITY( num_import, distributor.getNumImports() ); 
    TEST_EQUALITY( dim*distributor.getNumImports(), tgt_decomp_src.size() ); 
    for ( int i = 0; i < num_import; ++i )
    {
	TEST_EQUALITY( tgt_decomp_src[dim*i], 4.0+i );
	TEST_EQUALITY( tgt_decomp_src[dim*i+1], 2.0*inverse_rank );
    }

    Teuchos::Array<double> src_data(num_src_points);
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_data[i] = i*inverse_rank;
    }
    Teuchos::Array<double> tgt_data( distributor.getNumImports() );
    Teuchos::ArrayView<const double> src_view = src_data();
    distributor.distribute( src_view, tgt_data() );
    for ( int i = 0; i < num_import; ++i )
    {
	TEST_EQUALITY( tgt_data[i], (4.0+i)*rank );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CenterDistributor, dim_3_test )
{
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    int rank = comm->getRank();
    int size = comm->getSize();
    int inverse_rank = size - rank - 1;

    int dim = 3;
    int num_src_points = 10;
    int num_src_coords = dim*num_src_points;

    Teuchos::Array<double> src_coords(num_src_coords);
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_coords[dim*i] = 1.0*i;
	src_coords[dim*i+1] = 2.0*rank;
	src_coords[dim*i+2] = 2.0*rank;
    }

    int num_tgt_points = 2;
    int num_tgt_coords = dim*num_tgt_points;
    Teuchos::Array<double> tgt_coords( num_tgt_coords );
    tgt_coords[0] = 4.9;
    tgt_coords[1] = 2.0*inverse_rank;
    tgt_coords[2] = 2.0*inverse_rank;
    tgt_coords[3] = 11.4;
    tgt_coords[4] = 2.0*inverse_rank;
    tgt_coords[5] = 2.0*inverse_rank;

    double radius = 1.5;

    Teuchos::Array<double> tgt_decomp_src;

    DataTransferKit::CenterDistributor<3> distributor( 
	comm, src_coords(), tgt_coords(), radius, tgt_decomp_src );

    int num_import = 6;
    TEST_EQUALITY( num_import, distributor.getNumImports() ); 
    TEST_EQUALITY( dim*distributor.getNumImports(), tgt_decomp_src.size() );
    for ( int i = 0; i < num_import; ++i )
    {
	TEST_EQUALITY( tgt_decomp_src[dim*i], 4.0+i );
	TEST_EQUALITY( tgt_decomp_src[dim*i+1], 2.0*inverse_rank );
	TEST_EQUALITY( tgt_decomp_src[dim*i+2], 2.0*inverse_rank );
    }

    Teuchos::Array<double> src_data(num_src_points);
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_data[i] = i*inverse_rank;
    }
    Teuchos::Array<double> tgt_data( distributor.getNumImports() );
    Teuchos::ArrayView<const double> src_view = src_data();
    distributor.distribute( src_view, tgt_data() );
    for ( int i = 0; i < num_import; ++i )
    {
	TEST_EQUALITY( tgt_data[i], (4.0+i)*rank );
    }
}

//---------------------------------------------------------------------------//
// end tstCenterDistributor.cpp
//---------------------------------------------------------------------------//
