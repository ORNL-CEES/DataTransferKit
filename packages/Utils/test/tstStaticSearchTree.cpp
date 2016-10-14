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
 * \file   tstStaticSearchTree.cpp
 * \author Stuart R. Slattery
 * \brief  Static search tree unit tests.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <DTK_SearchTreeFactory.hpp>
#include <DTK_StaticSearchTree.hpp>

#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

// Get the default communicator.
Teuchos::RCP<const Teuchos::Comm<int>> getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<int>::getComm();
#else
    return Teuchos::rcp( new Teuchos::SerialComm<int>() );
#endif
}

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( NanoflannTree, dim_1_test )
{
    int dim = 1;
    int num_points = 10;
    int num_coords = dim * num_points;

    Teuchos::Array<double> coords( num_coords );
    for ( int i = 0; i < num_points; ++i )
    {
        coords[i] = 1.0 * i;
    }

    int max_leaf_size = 3;
    DataTransferKit::NanoflannTree<1> tree( coords(), max_leaf_size );

    Teuchos::Array<double> p1( dim );
    p1[0] = 4.9;
    Teuchos::Array<double> p2( dim );
    p2[0] = 11.4;

    int num_neighbors = 1;
    Teuchos::Array<unsigned> nnearest = tree.nnSearch( p1(), num_neighbors );
    TEST_EQUALITY( num_neighbors, nnearest.size() );
    TEST_EQUALITY( 5, nnearest[0] );

    nnearest = tree.nnSearch( p2(), num_neighbors );
    TEST_EQUALITY( num_neighbors, nnearest.size() );
    TEST_EQUALITY( 9, nnearest[0] );

    double radius = 1.1;
    nnearest = tree.radiusSearch( p1(), radius );
    TEST_EQUALITY( 3, nnearest.size() );
    TEST_EQUALITY( 5, nnearest[0] )
    TEST_EQUALITY( 4, nnearest[1] )
    TEST_EQUALITY( 6, nnearest[2] )

    nnearest = tree.radiusSearch( p2(), radius );
    TEST_EQUALITY( 0, nnearest.size() );

    radius = 2.41;
    nnearest = tree.radiusSearch( p2(), radius );
    TEST_EQUALITY( 1, nnearest.size() );
    TEST_EQUALITY( 9, nnearest[0] )
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( NanoflannTree, dim_1_factory_test )
{
    int dim = 1;
    int num_points = 10;
    int num_coords = dim * num_points;

    Teuchos::Array<double> coords( num_coords );
    for ( int i = 0; i < num_points; ++i )
    {
        coords[i] = 1.0 * i;
    }

    int max_leaf_size = 3;
    Teuchos::RCP<DataTransferKit::StaticSearchTree> tree =
        DataTransferKit::SearchTreeFactory::createStaticTree( 1, coords(),
                                                              max_leaf_size );

    Teuchos::Array<double> p1( dim );
    p1[0] = 4.9;
    Teuchos::Array<double> p2( dim );
    p2[0] = 11.4;

    int num_neighbors = 1;
    Teuchos::Array<unsigned> nnearest = tree->nnSearch( p1(), num_neighbors );
    TEST_EQUALITY( num_neighbors, nnearest.size() );
    TEST_EQUALITY( 5, nnearest[0] );

    nnearest = tree->nnSearch( p2(), num_neighbors );
    TEST_EQUALITY( num_neighbors, nnearest.size() );
    TEST_EQUALITY( 9, nnearest[0] );

    double radius = 1.1;
    nnearest = tree->radiusSearch( p1(), radius );
    TEST_EQUALITY( 3, nnearest.size() );
    TEST_EQUALITY( 5, nnearest[0] )
    TEST_EQUALITY( 4, nnearest[1] )
    TEST_EQUALITY( 6, nnearest[2] )

    nnearest = tree->radiusSearch( p2(), radius );
    TEST_EQUALITY( 0, nnearest.size() );

    radius = 2.41;
    nnearest = tree->radiusSearch( p2(), radius );
    TEST_EQUALITY( 1, nnearest.size() );
    TEST_EQUALITY( 9, nnearest[0] )
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( NanoflannTree, dim_2_test )
{
    int dim = 2;
    int num_points = 10;
    int num_coords = dim * num_points;

    Teuchos::Array<double> coords( num_coords );
    for ( int i = 0; i < num_points; ++i )
    {
        coords[dim * i] = 1.0 * i;
        coords[dim * i + 1] = 1.0;
    }

    int max_leaf_size = 2;
    DataTransferKit::NanoflannTree<2> tree( coords(), max_leaf_size );

    Teuchos::Array<double> p1( dim );
    p1[0] = 4.9;
    p1[1] = 1.0;
    Teuchos::Array<double> p2( dim );
    p2[0] = 11.4;
    p2[1] = 1.0;

    int num_neighbors = 1;
    Teuchos::Array<unsigned> nnearest = tree.nnSearch( p1(), num_neighbors );
    TEST_EQUALITY( num_neighbors, nnearest.size() );
    TEST_EQUALITY( 5, nnearest[0] );

    nnearest = tree.nnSearch( p2(), num_neighbors );
    TEST_EQUALITY( num_neighbors, nnearest.size() );
    TEST_EQUALITY( 9, nnearest[0] );

    double radius = 1.1;
    nnearest = tree.radiusSearch( p1(), radius );
    TEST_EQUALITY( 3, nnearest.size() );
    TEST_EQUALITY( 5, nnearest[0] )
    TEST_EQUALITY( 4, nnearest[1] )
    TEST_EQUALITY( 6, nnearest[2] )

    nnearest = tree.radiusSearch( p2(), radius );
    TEST_EQUALITY( 0, nnearest.size() );

    radius = 2.41;
    nnearest = tree.radiusSearch( p2(), radius );
    TEST_EQUALITY( 1, nnearest.size() );
    TEST_EQUALITY( 9, nnearest[0] )
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( NanoflannTree, dim_2_factory_test )
{
    int dim = 2;
    int num_points = 10;
    int num_coords = dim * num_points;

    Teuchos::Array<double> coords( num_coords );
    for ( int i = 0; i < num_points; ++i )
    {
        coords[dim * i] = 1.0 * i;
        coords[dim * i + 1] = 1.0;
    }

    int max_leaf_size = 2;
    Teuchos::RCP<DataTransferKit::StaticSearchTree> tree =
        DataTransferKit::SearchTreeFactory::createStaticTree( 2, coords(),
                                                              max_leaf_size );

    Teuchos::Array<double> p1( dim );
    p1[0] = 4.9;
    p1[1] = 1.0;
    Teuchos::Array<double> p2( dim );
    p2[0] = 11.4;
    p2[1] = 1.0;

    int num_neighbors = 1;
    Teuchos::Array<unsigned> nnearest = tree->nnSearch( p1(), num_neighbors );
    TEST_EQUALITY( num_neighbors, nnearest.size() );
    TEST_EQUALITY( 5, nnearest[0] );

    nnearest = tree->nnSearch( p2(), num_neighbors );
    TEST_EQUALITY( num_neighbors, nnearest.size() );
    TEST_EQUALITY( 9, nnearest[0] );

    double radius = 1.1;
    nnearest = tree->radiusSearch( p1(), radius );
    TEST_EQUALITY( 3, nnearest.size() );
    TEST_EQUALITY( 5, nnearest[0] )
    TEST_EQUALITY( 4, nnearest[1] )
    TEST_EQUALITY( 6, nnearest[2] )

    nnearest = tree->radiusSearch( p2(), radius );
    TEST_EQUALITY( 0, nnearest.size() );

    radius = 2.41;
    nnearest = tree->radiusSearch( p2(), radius );
    TEST_EQUALITY( 1, nnearest.size() );
    TEST_EQUALITY( 9, nnearest[0] )
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( NanoflannTree, dim_3_test )
{
    int dim = 3;
    int num_points = 10;
    int num_coords = dim * num_points;

    Teuchos::Array<double> coords( num_coords );
    for ( int i = 0; i < num_points; ++i )
    {
        coords[dim * i] = 1.0 * i;
        coords[dim * i + 1] = 1.0;
        coords[dim * i + 2] = 1.0;
    }

    int max_leaf_size = 3;
    DataTransferKit::NanoflannTree<3> tree( coords(), max_leaf_size );

    Teuchos::Array<double> p1( dim );
    p1[0] = 4.9;
    p1[1] = 1.0;
    p1[2] = 1.0;
    Teuchos::Array<double> p2( dim );
    p2[0] = 11.4;
    p2[1] = 1.0;
    p2[2] = 1.0;

    int num_neighbors = 1;
    Teuchos::Array<unsigned> nnearest = tree.nnSearch( p1(), num_neighbors );
    TEST_EQUALITY( num_neighbors, nnearest.size() );
    TEST_EQUALITY( 5, nnearest[0] );

    nnearest = tree.nnSearch( p2(), num_neighbors );
    TEST_EQUALITY( num_neighbors, nnearest.size() );
    TEST_EQUALITY( 9, nnearest[0] );

    double radius = 1.1;
    nnearest = tree.radiusSearch( p1(), radius );
    TEST_EQUALITY( 3, nnearest.size() );
    TEST_EQUALITY( 5, nnearest[0] )
    TEST_EQUALITY( 4, nnearest[1] )
    TEST_EQUALITY( 6, nnearest[2] )

    nnearest = tree.radiusSearch( p2(), radius );
    TEST_EQUALITY( 0, nnearest.size() );

    radius = 2.41;
    nnearest = tree.radiusSearch( p2(), radius );
    TEST_EQUALITY( 1, nnearest.size() );
    TEST_EQUALITY( 9, nnearest[0] )
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( NanoflannTree, dim_3_factory_test )
{
    int dim = 3;
    int num_points = 10;
    int num_coords = dim * num_points;

    Teuchos::Array<double> coords( num_coords );
    for ( int i = 0; i < num_points; ++i )
    {
        coords[dim * i] = 1.0 * i;
        coords[dim * i + 1] = 1.0;
        coords[dim * i + 2] = 1.0;
    }

    int max_leaf_size = 3;
    Teuchos::RCP<DataTransferKit::StaticSearchTree> tree =
        DataTransferKit::SearchTreeFactory::createStaticTree( 3, coords(),
                                                              max_leaf_size );

    Teuchos::Array<double> p1( dim );
    p1[0] = 4.9;
    p1[1] = 1.0;
    p1[2] = 1.0;
    Teuchos::Array<double> p2( dim );
    p2[0] = 11.4;
    p2[1] = 1.0;
    p2[2] = 1.0;

    int num_neighbors = 1;
    Teuchos::Array<unsigned> nnearest = tree->nnSearch( p1(), num_neighbors );
    TEST_EQUALITY( num_neighbors, nnearest.size() );
    TEST_EQUALITY( 5, nnearest[0] );

    nnearest = tree->nnSearch( p2(), num_neighbors );
    TEST_EQUALITY( num_neighbors, nnearest.size() );
    TEST_EQUALITY( 9, nnearest[0] );

    double radius = 1.1;
    nnearest = tree->radiusSearch( p1(), radius );
    TEST_EQUALITY( 3, nnearest.size() );
    TEST_EQUALITY( 5, nnearest[0] )
    TEST_EQUALITY( 4, nnearest[1] )
    TEST_EQUALITY( 6, nnearest[2] )

    nnearest = tree->radiusSearch( p2(), radius );
    TEST_EQUALITY( 0, nnearest.size() );

    radius = 2.41;
    nnearest = tree->radiusSearch( p2(), radius );
    TEST_EQUALITY( 1, nnearest.size() );
    TEST_EQUALITY( 9, nnearest[0] )
}

//---------------------------------------------------------------------------//
// end tstStaticSearchTree.cpp
//---------------------------------------------------------------------------//
