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
 * \file DTK_Box.hpp
 * \author Stuart R. Slattery
 * \brief Bounding box declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BOX_HPP
#define DTK_BOX_HPP

#include "DTK_BoundingBox.hpp"
#include "DTK_GeometryTraits.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_SerializationTraits.hpp>
#include <Teuchos_Tuple.hpp>

#include <iostream>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class Box
 * \brief Axis-aligned Cartesian box container.
 *
 * All three dimensions are explictly represented in this bounding box. This
 * is different from a bounding box in that it must always be finite and of a
 * fixed 3 dimensions.
 */
//---------------------------------------------------------------------------//
class Box
{

  public:
    // Default constructor.
    Box();

    // Constructor.
    Box( const double x_min, const double y_min, const double z_min,
         const double x_max, const double y_max, const double z_max );

    // Tuple constructor.
    Box( const Teuchos::Tuple<double, 6> &bounds );

    // Destructor.
    ~Box();

    // Determine if a point is in the box within a specified tolerance
    bool pointInBox( const Teuchos::Array<double> &coords,
                     const double tolerance ) const;

    // Get the boundaries of the box.
    Teuchos::Tuple<double, 6> getBounds() const
    {
        return Teuchos::tuple( d_x_min, d_y_min, d_z_min, d_x_max, d_y_max,
                               d_z_max );
    }

    // Compute the volume of the box.
    double volume() const;

    // Get the bounding box around the box.
    BoundingBox boundingBox() const;

    // Get the centroid of the box.
    Teuchos::Array<double> centroid() const;

  private:
    // X min.
    double d_x_min;

    // Y min.
    double d_y_min;

    // Z min.
    double d_z_min;

    // X max.
    double d_x_max;

    // Y max.
    double d_y_max;

    // Z max.
    double d_z_max;
};

//! overload for printing box
std::ostream &operator<<( std::ostream &os, const DataTransferKit::Box &b );

//---------------------------------------------------------------------------//
// GeometryTraits Specialization.
//---------------------------------------------------------------------------//
template <>
class GeometryTraits<Box>
{
  public:
    typedef Box geometry_type;

    static inline int dim( const Box & /*box*/ ) { return 3; }

    static inline double measure( const Box &box ) { return box.volume(); }

    static inline bool pointInGeometry( const Box &box,
                                        const Teuchos::Array<double> &coords,
                                        const double tolerance )
    {
        return box.pointInBox( coords, tolerance );
    }

    static inline BoundingBox boundingBox( const Box &box )
    {
        return box.boundingBox();
    }

    static inline Teuchos::Array<double> centroid( const Box &box )
    {
        return box.centroid();
    }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Serialization Traits Specialization.
//---------------------------------------------------------------------------//

namespace Teuchos
{

template <typename Ordinal>
class SerializationTraits<Ordinal, DataTransferKit::Box>
    : public DirectSerializationTraits<Ordinal, DataTransferKit::Box>
{
};

} // end namespace Teuchos

//---------------------------------------------------------------------------//

#endif // end DTK_BOX_HPP

//---------------------------------------------------------------------------//
// end DTK_Box.hpp
//---------------------------------------------------------------------------//
