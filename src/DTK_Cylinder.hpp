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
 * \file DTK_Cylinder.hpp
 * \author Stuart R. Slattery
 * \brief cylinder declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CYLINDER_HPP
#define DTK_CYLINDER_HPP

#include "DTK_BoundingBox.hpp"
#include "DTK_GeometryTraits.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_SerializationTraits.hpp>

#include <iostream>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class Cylinder
 * \brief Z-axis-aligned Cartesian cylinder container
 *
 * All three dimensions are explictly represented in this cylinder.
 */
//---------------------------------------------------------------------------//
class Cylinder
{

  public:

    // Default constructor.
    Cylinder();

    // Constructor.
    Cylinder( const double length, const double radius,
	      const double centroid_x, const double centroid_y,
	      const double centroid_z );

    // Destructor.
    ~Cylinder();

    // Determine if a point is in the cylinder within a specified tolerance.
    bool pointInCylinder( const Teuchos::Array<double>& coords,
			  const double tolerance ) const;

    //! Get the length of the cylinder.
    double length() const
    { return d_length; }

    //! Get the radius of the cylinder.
    double radius() const
    { return d_radius; }

    // Get the centroid of the cylinder.
    Teuchos::Array<double> centroid() const;

    // Compute the volume of the cylinder.
    double volume() const;

    // Compute the bounding box around the cylinder.
    BoundingBox boundingBox() const;

  private:
    
    // Length.
    double d_length;

    // Radius.
    double d_radius;

    // Centroid X-coordinate
    double d_centroid_x;

    // Centroid Y-coordinate
    double d_centroid_y;

    // Centroid Z-coordinate
    double d_centroid_z;
};

//! overload for printing cylinder
std::ostream& operator<< (std::ostream& os,const DataTransferKit::Cylinder& c); 

//---------------------------------------------------------------------------//
// GeometryTraits Specialization.
//---------------------------------------------------------------------------//
template<>
class GeometryTraits<Cylinder>
{
  public:

    typedef Cylinder geometry_type;

    static inline int dim( const Cylinder& cylinder )
    { return 3; }

    static inline double measure( const Cylinder& cylinder )
    { return cylinder.volume(); }

    static inline bool pointInGeometry( const Cylinder& cylinder,
					const Teuchos::Array<double>& coords,
					const double tolerance )
    { return cylinder.pointInCylinder( coords, tolerance ); }

    static inline BoundingBox boundingBox( const Cylinder& cylinder )
    { return cylinder.boundingBox(); }

    static inline Teuchos::Array<double> centroid( const Cylinder& cylinder )
    { return cylinder.centroid(); }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Serialization Traits Specialization.
//---------------------------------------------------------------------------//
namespace Teuchos
{

template<typename Ordinal>
class SerializationTraits<Ordinal, DataTransferKit::Cylinder>
    : public DirectSerializationTraits<Ordinal, DataTransferKit::Cylinder>
{};

} // end namespace Teuchos

//---------------------------------------------------------------------------//

#endif // end DTK_CYLINDER_HPP

//---------------------------------------------------------------------------//
// end DTK_Cylinder.hpp
//---------------------------------------------------------------------------//

