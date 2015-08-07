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
 * \file   tstCloudDomain.cpp
 * \author Stuart R. Slattery
 * \brief  CloudDomain unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <DTK_CloudDomain.hpp>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_Array.hpp"
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
TEUCHOS_UNIT_TEST( CloudDomain, dim_1_test )
{
    unsigned dim = 1;

    Teuchos::Array<double> domain_bounds(dim*2);
    domain_bounds[0] = 1.5;
    domain_bounds[1] = 1.7;

    DataTransferKit::CloudDomain<1> domain( domain_bounds.getRawPtr() );

    Teuchos::ArrayView<const double> bounds_view = domain.bounds();
    TEST_EQUALITY( 2*dim, bounds_view.size() );
    TEST_EQUALITY( domain_bounds[0], bounds_view[0] );
    TEST_EQUALITY( domain_bounds[1], bounds_view[1] );

    Teuchos::Array<double> center = domain.center();
    TEST_EQUALITY( dim, center.size() );
    TEST_EQUALITY( 1.6, center[0] );

    Teuchos::Array<double> p1( 1, 0.3 );
    Teuchos::Array<double> p2( 1, 1.55 );
    TEST_ASSERT( !domain.pointInDomain(p1()) );
    TEST_ASSERT( domain.pointInDomain(p2()) );

    double radius = 1.4;
    domain.expand( radius );
    TEST_EQUALITY( domain_bounds[0]-radius, bounds_view[0] );
    TEST_EQUALITY( domain_bounds[1]+radius, bounds_view[1] );
    TEST_EQUALITY( 1.6, center[0] );
    TEST_ASSERT( domain.pointInDomain(p1()) );
    TEST_ASSERT( domain.pointInDomain(p2()) );
}


//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CloudDomain, dim_2_test )
{
    int dim = 2;

    Teuchos::Array<double> domain_bounds(dim*2);
    domain_bounds[0] = 1.5;
    domain_bounds[1] = 1.7;
    domain_bounds[2] = 2.5;
    domain_bounds[3] = 2.7;

    DataTransferKit::CloudDomain<2> domain( domain_bounds.getRawPtr() );

    Teuchos::ArrayView<const double> bounds_view = domain.bounds();
    TEST_EQUALITY( 2*dim, bounds_view.size() );
    TEST_EQUALITY( domain_bounds[0], bounds_view[0] );
    TEST_EQUALITY( domain_bounds[1], bounds_view[1] );
    TEST_EQUALITY( domain_bounds[2], bounds_view[2] );
    TEST_EQUALITY( domain_bounds[3], bounds_view[3] );

    Teuchos::Array<double> center = domain.center();
    TEST_EQUALITY( dim, center.size() );
    TEST_EQUALITY( 1.6, center[0] );
    TEST_EQUALITY( 2.6, center[1] );

    Teuchos::Array<double> p1( 2 );
    p1[0] = 0.3;
    p1[1] = 1.3;
    Teuchos::Array<double> p2( 2 );
    p2[0] = 1.55;
    p2[1] = 2.55;
    TEST_ASSERT( !domain.pointInDomain(p1()) );
    TEST_ASSERT( domain.pointInDomain(p2()) );

    double radius = 1.4;
    domain.expand( radius );
    TEST_EQUALITY( domain_bounds[0]-radius, bounds_view[0] );
    TEST_EQUALITY( domain_bounds[1]+radius, bounds_view[1] );
    TEST_EQUALITY( domain_bounds[2]-radius, bounds_view[2] );
    TEST_EQUALITY( domain_bounds[3]+radius, bounds_view[3] );
    TEST_EQUALITY( 1.6, center[0] );
    TEST_EQUALITY( 2.6, center[1] );
    TEST_ASSERT( domain.pointInDomain(p1()) );
    TEST_ASSERT( domain.pointInDomain(p2()) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CloudDomain, dim_3_test )
{
    int dim = 3;


    Teuchos::Array<double> domain_bounds(dim*2);
    domain_bounds[0] = 1.5;
    domain_bounds[1] = 1.7;
    domain_bounds[2] = 2.5;
    domain_bounds[3] = 2.7;
    domain_bounds[4] = 3.5;
    domain_bounds[5] = 3.7;

    DataTransferKit::CloudDomain<3> domain( domain_bounds.getRawPtr() );

    Teuchos::ArrayView<const double> bounds_view = domain.bounds();
    TEST_EQUALITY( 2*dim, bounds_view.size() );
    TEST_EQUALITY( domain_bounds[0], bounds_view[0] );
    TEST_EQUALITY( domain_bounds[1], bounds_view[1] );
    TEST_EQUALITY( domain_bounds[2], bounds_view[2] );
    TEST_EQUALITY( domain_bounds[3], bounds_view[3] );
    TEST_EQUALITY( domain_bounds[4], bounds_view[4] );
    TEST_EQUALITY( domain_bounds[5], bounds_view[5] );

    Teuchos::Array<double> center = domain.center();
    TEST_EQUALITY( dim, center.size() );
    TEST_EQUALITY( 1.6, center[0] );
    TEST_EQUALITY( 2.6, center[1] );
    TEST_EQUALITY( 3.6, center[2] );

    Teuchos::Array<double> p1( 3 );
    p1[0] = 0.3;
    p1[1] = 1.3;
    p1[2] = 2.3;
    Teuchos::Array<double> p2( 3 );
    p2[0] = 1.55;
    p2[1] = 2.55;
    p2[2] = 3.55;
    TEST_ASSERT( !domain.pointInDomain(p1()) );
    TEST_ASSERT( domain.pointInDomain(p2()) );

    double radius = 1.4;
    domain.expand( radius );
    TEST_EQUALITY( domain_bounds[0]-radius, bounds_view[0] );
    TEST_EQUALITY( domain_bounds[1]+radius, bounds_view[1] );
    TEST_EQUALITY( domain_bounds[2]-radius, bounds_view[2] );
    TEST_EQUALITY( domain_bounds[3]+radius, bounds_view[3] );
    TEST_EQUALITY( domain_bounds[4]-radius, bounds_view[4] );
    TEST_EQUALITY( domain_bounds[5]+radius, bounds_view[5] );
    TEST_EQUALITY( 1.6, center[0] );
    TEST_EQUALITY( 2.6, center[1] );
    TEST_EQUALITY( 3.6, center[2] );
    TEST_ASSERT( domain.pointInDomain(p1()) );
    TEST_ASSERT( domain.pointInDomain(p2()) );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( CloudDomain, dim_3_parallel_test )
{
    int dim = 3;

    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

    Teuchos::Array<double> domain_bounds(dim*2);
    domain_bounds[0] = 1.5;
    domain_bounds[1] = 1.7;
    domain_bounds[2] = 2.5;
    domain_bounds[3] = 2.7;
    domain_bounds[4] = 3.5;
    domain_bounds[5] = 3.7;

    DataTransferKit::CloudDomain<3> domain;
    if ( 0 == comm->getRank() )
    {
	domain = DataTransferKit::CloudDomain<3>( domain_bounds.getRawPtr() );
    }

    Teuchos::broadcast( *comm, 0, 
			Teuchos::Ptr<DataTransferKit::CloudDomain<3> >(&domain) );

    Teuchos::ArrayView<const double> bounds_view = domain.bounds();
    TEST_EQUALITY( 2*dim, bounds_view.size() );
    TEST_EQUALITY( domain_bounds[0], bounds_view[0] );
    TEST_EQUALITY( domain_bounds[1], bounds_view[1] );
    TEST_EQUALITY( domain_bounds[2], bounds_view[2] );
    TEST_EQUALITY( domain_bounds[3], bounds_view[3] );
    TEST_EQUALITY( domain_bounds[4], bounds_view[4] );
    TEST_EQUALITY( domain_bounds[5], bounds_view[5] );

    Teuchos::Array<double> center = domain.center();
    TEST_EQUALITY( dim, center.size() );
    TEST_EQUALITY( 1.6, center[0] );
    TEST_EQUALITY( 2.6, center[1] );
    TEST_EQUALITY( 3.6, center[2] );

    Teuchos::Array<double> p1( 3 );
    p1[0] = 0.3;
    p1[1] = 1.3;
    p1[2] = 2.3;
    Teuchos::Array<double> p2( 3 );
    p2[0] = 1.55;
    p2[1] = 2.55;
    p2[2] = 3.55;
    TEST_ASSERT( !domain.pointInDomain(p1()) );
    TEST_ASSERT( domain.pointInDomain(p2()) );
}

//---------------------------------------------------------------------------//
// end tstCloudDomain.cpp
//---------------------------------------------------------------------------//
