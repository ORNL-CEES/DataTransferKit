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
 * \file   tstSplineInterpolationPairing.cpp
 * \author Stuart R. Slattery
 * \brief  Interpolation pairing.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <DTK_SplineInterpolationPairing.hpp>

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
// TEST EPSILON
//---------------------------------------------------------------------------//

const double epsilon = 1.0e-14;

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( SplineInterpolationPairing, radius_dim_1_test )
{
    int dim = 1;
    int num_src_points = 10;
    int num_src_coords = dim*num_src_points;

    Teuchos::Array<double> src_coords(num_src_coords);
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_coords[dim*i] = 1.0*i;
    }

    int num_tgt_points = 2;
    int num_tgt_coords = dim*num_tgt_points;
    Teuchos::Array<double> tgt_coords( num_tgt_coords );
    tgt_coords[0] = 4.9;
    tgt_coords[1] = 10.0;

    double radius = 1.1;

    DataTransferKit::SplineInterpolationPairing<1> pairing( 
	src_coords(), tgt_coords(), false, 0, radius );
    
    Teuchos::ArrayView<const unsigned> view = pairing.childCenterIds( 0 );
    TEST_EQUALITY( 3, view.size() );
    TEST_EQUALITY( 5, view[0] )
    TEST_EQUALITY( 4, view[1] )
    TEST_EQUALITY( 6, view[2] )

    view = pairing.childCenterIds( 1 );
    TEST_EQUALITY( 1, view.size() );
    TEST_EQUALITY( 9, view[0] );

    Teuchos::ArrayRCP<DataTransferKit::EntityId> children_per_parent = 
	pairing.childrenPerParent();
    TEST_EQUALITY( children_per_parent[0], 3 );
    TEST_EQUALITY( children_per_parent[1], 1 );

    double parent_radius = pairing.parentSupportRadius( 0 );
    TEST_EQUALITY( parent_radius, radius );

    parent_radius = pairing.parentSupportRadius( 1 );
    TEST_EQUALITY( parent_radius, radius );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( SplineInterpolationPairing, radius_dim_2_test )
{
    int dim = 2;
    int num_src_points = 10;
    int num_src_coords = dim*num_src_points;

    Teuchos::Array<double> src_coords(num_src_coords);
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_coords[dim*i] = 1.0*i;
	src_coords[dim*i+1] = 1.0;
    }

    int num_tgt_points = 2;
    int num_tgt_coords = dim*num_tgt_points;
    Teuchos::Array<double> tgt_coords( num_tgt_coords );
    tgt_coords[0] = 4.9;
    tgt_coords[1] = 1.0;
    tgt_coords[2] = 10.0;
    tgt_coords[3] = 1.0;

    double radius = 1.1;

    DataTransferKit::SplineInterpolationPairing<2> pairing( 
	src_coords(), tgt_coords(), false, 0, radius );
    
    Teuchos::ArrayView<const unsigned> view = pairing.childCenterIds( 0 );
    TEST_EQUALITY( 3, view.size() );
    TEST_EQUALITY( 5, view[0] )
    TEST_EQUALITY( 4, view[1] )
    TEST_EQUALITY( 6, view[2] )

    view = pairing.childCenterIds( 1 );
    TEST_EQUALITY( 1, view.size() );
    TEST_EQUALITY( 9, view[0] );

    Teuchos::ArrayRCP<DataTransferKit::EntityId> children_per_parent = 
	pairing.childrenPerParent();
    TEST_EQUALITY( children_per_parent[0], 3 );
    TEST_EQUALITY( children_per_parent[1], 1 );

    double parent_radius = pairing.parentSupportRadius( 0 );
    TEST_EQUALITY( parent_radius, radius );

    parent_radius = pairing.parentSupportRadius( 1 );
    TEST_EQUALITY( parent_radius, radius );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( SplineInterpolationPairing, radius_dim_3_test )
{
    int dim = 3;
    int num_src_points = 10;
    int num_src_coords = dim*num_src_points;

    Teuchos::Array<double> src_coords(num_src_coords);
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_coords[dim*i] = 1.0*i;
	src_coords[dim*i+1] = 1.0;
	src_coords[dim*i+2] = 1.0;
    }

    int num_tgt_points = 2;
    int num_tgt_coords = dim*num_tgt_points;
    Teuchos::Array<double> tgt_coords( num_tgt_coords );
    tgt_coords[0] = 4.9;
    tgt_coords[1] = 1.0;
    tgt_coords[2] = 1.0;
    tgt_coords[3] = 10.0;
    tgt_coords[4] = 1.0;
    tgt_coords[5] = 1.0;

    double radius = 1.1;

    DataTransferKit::SplineInterpolationPairing<3> pairing( 
	src_coords(), tgt_coords(), false, 0, radius );
    
    Teuchos::ArrayView<const unsigned> view = pairing.childCenterIds( 0 );
    TEST_EQUALITY( 3, view.size() );
    TEST_EQUALITY( 5, view[0] )
    TEST_EQUALITY( 4, view[1] )
    TEST_EQUALITY( 6, view[2] )

    view = pairing.childCenterIds( 1 );
    TEST_EQUALITY( 1, view.size() );
    TEST_EQUALITY( 9, view[0] );

    Teuchos::ArrayRCP<DataTransferKit::EntityId> children_per_parent = 
	pairing.childrenPerParent();
    TEST_EQUALITY( children_per_parent[0], 3 );
    TEST_EQUALITY( children_per_parent[1], 1 );

    double parent_radius = pairing.parentSupportRadius( 0 );
    TEST_EQUALITY( parent_radius, radius );

    parent_radius = pairing.parentSupportRadius( 1 );
    TEST_EQUALITY( parent_radius, radius );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( SplineInterpolationPairing, knn_dim_1_test )
{
    int dim = 1;
    int num_src_points = 10;
    int num_src_coords = dim*num_src_points;

    Teuchos::Array<double> src_coords(num_src_coords);
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_coords[dim*i] = 1.0*i;
    }

    int num_tgt_points = 2;
    int num_tgt_coords = dim*num_tgt_points;
    Teuchos::Array<double> tgt_coords( num_tgt_coords );
    tgt_coords[0] = 4.9;
    tgt_coords[1] = 10.0;

    unsigned knn = 3;

    DataTransferKit::SplineInterpolationPairing<1> pairing( 
	src_coords(), tgt_coords(), true, knn, 0.0 );
    
    Teuchos::ArrayView<const unsigned> view = pairing.childCenterIds( 0 );
    TEST_EQUALITY( knn, view.size() );
    TEST_EQUALITY( 5, view[0] );
    TEST_EQUALITY( 4, view[1] );
    TEST_EQUALITY( 6, view[2] );

    view = pairing.childCenterIds( 1 );
    TEST_EQUALITY( knn, view.size() );
    TEST_EQUALITY( 9, view[0] );
    TEST_EQUALITY( 8, view[1] );
    TEST_EQUALITY( 7, view[2] );

    Teuchos::ArrayRCP<DataTransferKit::EntityId> children_per_parent = 
	pairing.childrenPerParent();
    TEST_EQUALITY( children_per_parent[0], knn );
    TEST_EQUALITY( children_per_parent[1], knn );

    double radius = pairing.parentSupportRadius( 0 );
    TEST_FLOATING_EQUALITY( 1.01*1.1, radius, epsilon );

    radius = pairing.parentSupportRadius( 1 );
    TEST_FLOATING_EQUALITY( 1.01*3.0, radius, epsilon );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( SplineInterpolationPairing, knn_dim_2_test )
{
    int dim = 2;
    int num_src_points = 10;
    int num_src_coords = dim*num_src_points;

    Teuchos::Array<double> src_coords(num_src_coords);
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_coords[dim*i] = 1.0*i;
	src_coords[dim*i+1] = 1.0;
    }

    int num_tgt_points = 2;
    int num_tgt_coords = dim*num_tgt_points;
    Teuchos::Array<double> tgt_coords( num_tgt_coords );
    tgt_coords[0] = 4.9;
    tgt_coords[1] = 1.0;
    tgt_coords[2] = 10.0;
    tgt_coords[3] = 1.0;

    unsigned knn = 3;

    DataTransferKit::SplineInterpolationPairing<2> pairing( 
	src_coords(), tgt_coords(), true, knn, 0.0 );
    
    Teuchos::ArrayView<const unsigned> view = pairing.childCenterIds( 0 );
    TEST_EQUALITY( knn, view.size() );
    TEST_EQUALITY( 5, view[0] );
    TEST_EQUALITY( 4, view[1] );
    TEST_EQUALITY( 6, view[2] );

    view = pairing.childCenterIds( 1 );
    TEST_EQUALITY( knn, view.size() );
    TEST_EQUALITY( 9, view[0] );
    TEST_EQUALITY( 8, view[1] );
    TEST_EQUALITY( 7, view[2] );

    Teuchos::ArrayRCP<DataTransferKit::EntityId> children_per_parent = 
	pairing.childrenPerParent();
    TEST_EQUALITY( children_per_parent[0], knn );
    TEST_EQUALITY( children_per_parent[1], knn );

    double radius = pairing.parentSupportRadius( 0 );
    TEST_FLOATING_EQUALITY( 1.01*1.1, radius, epsilon );

    radius = pairing.parentSupportRadius( 1 );
    TEST_FLOATING_EQUALITY( 1.01*3.0, radius, epsilon );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( SplineInterpolationPairing, knn_dim_3_test )
{
    int dim = 3;
    int num_src_points = 10;
    int num_src_coords = dim*num_src_points;

    Teuchos::Array<double> src_coords(num_src_coords);
    for ( int i = 0; i < num_src_points; ++i )
    {
	src_coords[dim*i] = 1.0*i;
	src_coords[dim*i+1] = 1.0;
	src_coords[dim*i+2] = 1.0;
    }

    int num_tgt_points = 2;
    int num_tgt_coords = dim*num_tgt_points;
    Teuchos::Array<double> tgt_coords( num_tgt_coords );
    tgt_coords[0] = 4.9;
    tgt_coords[1] = 1.0;
    tgt_coords[2] = 1.0;
    tgt_coords[3] = 10.0;
    tgt_coords[4] = 1.0;
    tgt_coords[5] = 1.0;

    unsigned knn = 3;

    DataTransferKit::SplineInterpolationPairing<3> pairing( 
	src_coords(), tgt_coords(), true, knn, 0.0 );
    
    Teuchos::ArrayView<const unsigned> view = pairing.childCenterIds( 0 );
    TEST_EQUALITY( knn, view.size() );
    TEST_EQUALITY( 5, view[0] );
    TEST_EQUALITY( 4, view[1] );
    TEST_EQUALITY( 6, view[2] );

    view = pairing.childCenterIds( 1 );
    TEST_EQUALITY( knn, view.size() );
    TEST_EQUALITY( 9, view[0] );
    TEST_EQUALITY( 8, view[1] );
    TEST_EQUALITY( 7, view[2] );

    Teuchos::ArrayRCP<DataTransferKit::EntityId> children_per_parent = 
	pairing.childrenPerParent();
    TEST_EQUALITY( children_per_parent[0], knn );
    TEST_EQUALITY( children_per_parent[1], knn );

    double radius = pairing.parentSupportRadius( 0 );
    TEST_FLOATING_EQUALITY( 1.01*1.1, radius, epsilon );

    radius = pairing.parentSupportRadius( 1 );
    TEST_FLOATING_EQUALITY( 1.01*3.0, radius, epsilon );
}

//---------------------------------------------------------------------------//
// end tstSplineInterpolationPairing.cpp
//---------------------------------------------------------------------------//
