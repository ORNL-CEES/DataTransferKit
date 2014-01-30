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
  * \file DTK_Box.cpp
  * \author Stuart R. Slattery
  * \brief Bounding box definition.
  */
//---------------------------------------------------------------------------//

#include "DTK_Box.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
Box::Box()
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
 *
 * \param x_min Minimum x coordinate value in the box.
 *
 * \param y_min Minimum y coordinate value in the box.
 *
 * \param z_min Minimum z coordinate value in the box.
 *
 * \param x_max Maximum x coordinate value in the box.
 *
 * \param y_max Maximum y coordinate value in the box.
 *
 * \param z_max Maximum z coordinate value in the box.
 */
Box::Box( 
    const double x_min, const double y_min, const double z_min,
    const double x_max, const double y_max, const double z_max )
    : d_x_min( x_min )
    , d_y_min( y_min )
    , d_z_min( z_min )
    , d_x_max( x_max )
    , d_y_max( y_max )
    , d_z_max( z_max )
{
    DTK_REQUIRE( d_x_min <= d_x_max );
    DTK_REQUIRE( d_y_min <= d_y_max );
    DTK_REQUIRE( d_z_min <= d_z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tuple constructor.
 *
 * \param bounds Tuple containing {x_min, y_min, z_min, x_max, y_max, z_max}.
 */
Box::Box( const Teuchos::Tuple<double,6>& bounds )
    : d_x_min( bounds[0] )
    , d_y_min( bounds[1] )
    , d_z_min( bounds[2] )
    , d_x_max( bounds[3] )
    , d_y_max( bounds[4] )
    , d_z_max( bounds[5] )
{
    DTK_REQUIRE( d_x_min <= d_x_max );
    DTK_REQUIRE( d_y_min <= d_y_max );
    DTK_REQUIRE( d_z_min <= d_z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Box::~Box()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if a point is in the box within a specified tolerance.
 *
 * \param coords Cartesian coordinates to check for point inclusion. The
 * coordinates must have a dimension between 0 and 3.
 *
 * \param tolerance The geometric tolerance to check point-inclusion with.
 *
 * \return Return true if the point is in the box, false if not. A point on
 * the box boundary will return true.
 */
bool Box::pointInBox( const Teuchos::Array<double>& coords,
		      const double tolerance ) const
{
    DTK_REQUIRE( 3 == coords.size() );

    if ( coords[0] >= d_x_min - tolerance &&
	 coords[1] >= d_y_min - tolerance &&
	 coords[2] >= d_z_min - tolerance &&
	 coords[0] <= d_x_max + tolerance &&
	 coords[1] <= d_y_max + tolerance &&
	 coords[2] <= d_z_max + tolerance )
    {
	return true;
    }

    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the volume of the box.
 *
 * \return Return the volume of the box.
 */
double Box::volume() const
{
    return (d_x_max-d_x_min)*(d_y_max-d_y_min)*(d_z_max-d_z_min);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the bounding box around the box.
 *
 * \return The bounding box around the box.
 */
BoundingBox Box::boundingBox() const
{
    return BoundingBox( d_x_min, d_y_min, d_z_min, d_x_max, d_y_max, d_z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the centroid of the box.
 *
 * \return The centroid coordinates.
 */
Teuchos::Array<double> Box::centroid() const
{
    Teuchos::Array<double> center(3);
    center[0] = (d_x_max + d_x_min) / 2.0;
    center[1] = (d_y_max + d_y_min) / 2.0;
    center[2] = (d_z_max + d_z_min) / 2.0;
    return center;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Print the box description to an ostream.
 *
 * \return The ostream.
 */

//---------------------------------------------------------------------------//
std::ostream& operator<< (std::ostream& os,const DataTransferKit::Box& b)
{
  Teuchos::Tuple<double,6> bounds = b.getBounds();

  os << "Box: d_x_min=" << bounds[0] 
     << ",d_y_min=" << bounds[1]
     << ",d_z_min=" << bounds[2]
     << ",d_x_max=" << bounds[3]
     << ",d_y_max=" << bounds[4]
     << ",d_z_max=" << bounds[5];

  return os;
}

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_Box.cpp
//---------------------------------------------------------------------------//

