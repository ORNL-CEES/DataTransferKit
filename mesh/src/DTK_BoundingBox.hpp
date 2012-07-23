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
 * \file DTK_BoundingBox.hpp
 * \author Stuart R. Slattery
 * \brief Bounding box declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BOUNDINGBOX_HPP
#define DTK_BOUNDINGBOX_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_SerializationTraits.hpp>

namespace DataTransferKit
{

//---------------------------------------------------------------------------//
/*!
 * \class BoundingBox
 * \brief Axis-aligned Cartesian bounding box container
 *
 * All three dimensions are explictly represented in this bounding box,
 * however, from an algorithmic standpoint, this can be treated as a one or
 * two dimensional box as well by setting the unused dimension bounds to +/-
 * Teuchos::ScalarTraits<double>::rmax().
 */
//---------------------------------------------------------------------------//
class BoundingBox
{

  public:

    // Default constructor.
    BoundingBox();

    // Constructor.
    BoundingBox( const double x_min, const double y_min, const double z_min,
		 const double x_max, const double y_max, const double z_max );

    // Tuple constructor.
    BoundingBox( const Teuchos::Tuple<double,6>& bounds );

    // Destructor.
    ~BoundingBox();

    // Determine if a point is in the box.
    bool pointInBox( const Teuchos::Array<double>& coords ) const;

    // Get the boundaries of the box.
    Teuchos::Tuple<double,6> getBounds() const
    { return Teuchos::tuple( d_x_min, d_y_min, d_z_min, 
			     d_x_max, d_y_max, d_z_max ); }

    // Compute the volume of the box given its dimension.
    double volume( const int dim ) const;

    // Static function for box intersection.
    static bool intersectBoxes( const BoundingBox& box_A,
				const BoundingBox& box_B,
				BoundingBox& intersection );
    
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

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_BOUNDINGBOX_HPP

//---------------------------------------------------------------------------//
// end DTK_BoundingBox.hpp
//---------------------------------------------------------------------------//

