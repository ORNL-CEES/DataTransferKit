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
  * \file DTK_BoundingBox.cpp
  * \author Stuart R. Slattery
  * \brief Bounding box definition.
  */
//---------------------------------------------------------------------------//

#include "DTK_BoundingBox.hpp"
#include <DTK_Exception.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
BoundingBox::BoundingBox()
    : d_x_min( 0.0 )
    , d_y_min( 0.0 )
    , d_z_min( 0.0 )
    , d_x_max( 0.0 )
    , d_y_max( 0.0 )
    , d_z_max( 0.0 )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
BoundingBox::BoundingBox( 
    const double x_min, const double y_min, const double z_min,
    const double x_max, const double y_max, const double z_max )
    : d_x_min( x_min )
    , d_y_min( y_min )
    , d_z_min( z_min )
    , d_x_max( x_max )
    , d_y_max( y_max )
    , d_z_max( z_max )
{
    testPrecondition( d_x_min <= d_x_max );
    testPrecondition( d_y_min <= d_y_max );
    testPrecondition( d_z_min <= d_z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tuple constructor.
 */
BoundingBox::BoundingBox( const Teuchos::Tuple<double,6>& bounds )
    : d_x_min( bounds[0] )
    , d_y_min( bounds[1] )
    , d_z_min( bounds[2] )
    , d_x_max( bounds[3] )
    , d_y_max( bounds[4] )
    , d_z_max( bounds[5] )
{
    testPrecondition( d_x_min <= d_x_max );
    testPrecondition( d_y_min <= d_y_max );
    testPrecondition( d_z_min <= d_z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
BoundingBox::~BoundingBox()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if a point is in the box. A point on the box boundary will
 * return true.
 */
bool BoundingBox::pointInBox( const Teuchos::Array<double>& coords ) const
{
    if ( coords.size() == 1 )
    {
	if ( coords[0] >= d_x_min &&
	     coords[0] <= d_x_max )
	{
	    return true;
	}
    }

    else if ( coords.size() == 2 )
    {
	if ( coords[0] >= d_x_min &&
	     coords[1] >= d_y_min &&
	     coords[0] <= d_x_max &&
	     coords[1] <= d_y_max )
	{
	    return true;
	}
    }

    else if ( coords.size() == 3 )
    {
	if ( coords[0] >= d_x_min &&
	     coords[1] >= d_y_min &&
	     coords[2] >= d_z_min &&
	     coords[0] <= d_x_max &&
	     coords[1] <= d_y_max &&
	     coords[2] <= d_z_max )
	{
	    return true;
	}
    }

    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the volume of the bounding box given its dimension.
 */
double BoundingBox::volume( const int dim ) const
{
    if ( dim == 1 )
    { 
	return (d_x_max-d_x_min);
    }
    else if ( dim == 2 )
    { 
	return (d_x_max-d_x_min)*(d_y_max-d_y_min);
    }
    else if ( dim == 3 )
    { 
	return (d_x_max-d_x_min)*(d_y_max-d_y_min)*(d_z_max-d_z_min);
    }
    else
    {
	testPrecondition( dim < 4 );
	return 0;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Static function for box intersection. Return false if they do not
 * intersect. 
 */
bool BoundingBox::intersectBoxes( const BoundingBox& box_A,
				  const BoundingBox& box_B,
				  BoundingBox& intersection)
{
    Teuchos::Tuple<double,6> bounds_A = box_A.getBounds();
    Teuchos::Tuple<double,6> bounds_B = box_B.getBounds();

    double x_min, y_min, z_min, x_max, y_max, z_max;

    // Test for no overlap in X.
    if ( bounds_A[0] > bounds_B[3] || bounds_A[3] < bounds_B[0] )
    {
	return false;
    }
    // Test for no overlap in Y.
    if ( bounds_A[1] > bounds_B[4] || bounds_A[4] < bounds_B[1] )
    {
	return false;
    }
    // Test for no overlap in Z.
    if ( bounds_A[2] > bounds_B[5] || bounds_A[5] < bounds_B[2] )
    {
	return false;
    }

    // Get overlap in X.
    if ( bounds_A[0] > bounds_B[0] )
    {
	x_min = bounds_A[0];
    }
    else
    {
	x_min = bounds_B[0];
    }
    if ( bounds_A[3] > bounds_B[3] )
    {
	x_max = bounds_B[3];
    }
    else
    {
	x_max = bounds_A[3];
    }

    // Get overlap in Y.
    if ( bounds_A[1] > bounds_B[1] )
    {
	y_min = bounds_A[1];
    }
    else
    {
	y_min = bounds_B[1];
    }
    if ( bounds_A[4] > bounds_B[4] )
    {
	y_max = bounds_B[4];
    }
    else
    {
	y_max = bounds_A[4];
    }

    // Get overlap in Z.
    if ( bounds_A[2] > bounds_B[2] )
    {
	z_min = bounds_A[2];
    }
    else
    {
	z_min = bounds_B[2];
    }
    if ( bounds_A[5] > bounds_B[5] )
    {
	z_max = bounds_B[5];
    }
    else
    {
	z_max = bounds_A[5];
    }

    intersection = BoundingBox( x_min, y_min, z_min, x_max, y_max, z_max );
    return true;
}

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_BoundingBox.cpp
//---------------------------------------------------------------------------//

