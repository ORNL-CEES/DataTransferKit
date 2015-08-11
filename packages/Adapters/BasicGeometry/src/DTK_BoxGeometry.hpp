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
 * \file DTK_BoxGeometry.hpp
 * \author Stuart R. Slattery
 * \brief BoxGeometry declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BOXGEOMETRY_HPP
#define DTK_BOXGEOMETRY_HPP

#include "DTK_BasicGeometryEntity.hpp"

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ArrayView.hpp>

#include <iostream>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class BoxGeometry
 * \brief Axis-aligned Cartesian box container.
 *
 * All three dimensions are explictly represented in this bounding box. This
 * is different from a bounding box in that it must always be finite and of a
 * fixed 3 dimensions.
 */
//---------------------------------------------------------------------------//
class BoxGeometry : public BasicGeometryEntity
{
  public:

    // Default constructor.
    BoxGeometry();

    // Constructor.
    BoxGeometry( const EntityId global_id, const int owner_rank, const int block_id,
	 const double x_min, const double y_min, const double z_min,
	 const double x_max, const double y_max, const double z_max );

    // Tuple constructor.
    BoxGeometry( const EntityId global_id,
	 const int owner_rank, 
	 const int block_id,
	 const Teuchos::Tuple<double,6>& bounds );

    // Destructor.
    ~BoxGeometry();

    // Static function to check for box intersection but not perform it.
    static bool checkForIntersection( const BoxGeometry& box_A,
				      const BoxGeometry& box_B );

    // Static function for box intersection.
    static bool intersectBoxGeometryes( const BoxGeometry& box_A,
				const BoxGeometry& box_B,
				BoxGeometry& box_intersection );

    // Static function for box union
    static void uniteBoxGeometryes( const BoxGeometry& box_A,
			    const BoxGeometry& box_B,
			    BoxGeometry& box_union );

    // Compound assignment overload.
    BoxGeometry& operator+=(const BoxGeometry& rhs);

    // Return the entity measure.
    double measure() const override;

    // Compute the centroid of the entity.
    void centroid( const Teuchos::ArrayView<double>& centroid ) const override;

    // (Reverse Map) Map a point to the reference space of an entity. Return
    // the parameterized point.
    bool mapToReferenceFrame( 
	const Teuchos::ArrayView<const double>& point,
	const Teuchos::ArrayView<double>& reference_point ) const override;

    // Determine if a reference point is in the parameterized space of an
    // entity.
    bool checkPointInclusion( 
	const double tolerance,
	const Teuchos::ArrayView<const double>& reference_point ) const override;

    // (Forward Map) Map a reference point to the physical space of an entity.
    void mapToPhysicalFrame( 
	const Teuchos::ArrayView<const double>& reference_point,
	const Teuchos::ArrayView<double>& point ) const override;
};

//---------------------------------------------------------------------------//
//! Addition operator overload. Adding two boxes together is the same as
//! computing their union.
BoxGeometry operator+( const BoxGeometry& box_1, const BoxGeometry& box_2 );

//---------------------------------------------------------------------------//
//! Overload for printing box.
std::ostream& operator<< (std::ostream& os,const DataTransferKit::BoxGeometry& b); 

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_BOXGEOMETRY_HPP

//---------------------------------------------------------------------------//
// end DTK_BoxGeometry.hpp
//---------------------------------------------------------------------------//

