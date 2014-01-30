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
  * \file DTK_Cylinder.cpp
  * \author Stuart R. Slattery
  * \brief Bounding box definition.
  */
//---------------------------------------------------------------------------//

#include <cmath>

#include "DTK_Cylinder.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
Cylinder::Cylinder()
    : d_length( 0.0 )
    , d_radius( 0.0 )
    , d_centroid_x( 0.0 )
    , d_centroid_y( 0.0 )
    , d_centroid_z( 0.0 )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param length Length of cylinder along Z-axis.
 *
 * \param radius Radius of cylinder.
 *
 * \param centroid_x Centroid X-coordinate.
 *
 * \param centroid_y Centroid Y-coordinate.
 *
 * \param centroid_z Centroid Z-coordinate.
 */
Cylinder::Cylinder( const double length, const double radius,
		    const double centroid_x, const double centroid_y,
		    const double centroid_z )
    : d_length( length )
    , d_radius( radius )
    , d_centroid_x( centroid_x )
    , d_centroid_y( centroid_y )
    , d_centroid_z( centroid_z )
{
    DTK_REQUIRE( 0.0 <= d_length );
    DTK_REQUIRE( 0.0 <= d_radius );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Cylinder::~Cylinder()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if a point is in the cylinder. 
 *
 * \param coords Cartesian coordinates to check for point inclusion. The
 * coordinates must have a dimension of 3.
 *
 * \param tolerance The geometric tolerance to check point-inclusion with.
 *
 * \return Return true if the point is in the cylinder, false if not. A point
 * on the cylinder boundary or outside but within the tolerance will return
 * true.
 */
bool Cylinder::pointInCylinder( const Teuchos::Array<double>& coords,
				const double tolerance ) const
{
    DTK_REQUIRE( coords.size() == 3 );

    double distance = std::pow(
	(d_centroid_x - coords[0])*(d_centroid_x - coords[0]) +
	(d_centroid_y - coords[1])*(d_centroid_y - coords[1]),
	0.5 );

    if ( distance <= d_radius + tolerance &&
	 coords[2] >= d_centroid_z - d_length/2 - tolerance &&
	 coords[2] <= d_centroid_z + d_length/2 + tolerance )
    {
	return true;
    }

    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the centroid of the cylinder.
 *
 * \return The centroid coordinates.
 */
Teuchos::Array<double> Cylinder::centroid() const
{
    Teuchos::Array<double> coords(3);
    coords[0] = d_centroid_x;
    coords[1] = d_centroid_y;
    coords[2] = d_centroid_z;
    return coords;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the volume of the cylinder given its dimension.
 *
 * \param dim The dimension of the cylinder we want to compute the volume
 * for. We need this because the cylinder always stores all 3
 * dimensions. Lower dimension cylinderes are resolved with higher dimensions
 * set to +/- Teuchos::ScalarTraits<double>::rmax(). For dim = 1, only the x
 * dimension is used. For dim = 2, the x and y dimensions are used. For dim =
 * 3, the x, y, and z dimensions are used.
 *
 * \return Return the volume of the cylinder.
 */
double Cylinder::volume() const
{
    double zero = 0.0;
    double pi = 2.0 * std::acos( zero );
    return pi * d_radius * d_radius * d_length;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the bounding box around the cylinder.
 *
 * \return The bounding box around the cylinder.
 */
BoundingBox Cylinder::boundingBox() const
{
    return BoundingBox( d_centroid_x - d_radius,
			d_centroid_y - d_radius,
			d_centroid_z - d_length/2,
			d_centroid_x + d_radius,
			d_centroid_y + d_radius,
			d_centroid_z + d_length/2 );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Print the cylinder description to an ostream.
 *
 * \return The ostream.
 */

//---------------------------------------------------------------------------//
std::ostream& operator<< (std::ostream& os,const DataTransferKit::Cylinder& c)
{
  os << "Cylinder: length=" << c.length() << ",radius=" << c.radius()
     << ", centroid=(" << c.centroid()[0] << "," << c.centroid()[1]
     << "," << c.centroid()[2] << ")";

  return os;
}

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_Cylinder.cpp
//---------------------------------------------------------------------------//

