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

#include "DTK_BasicGeometryEntity.hpp"

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ArrayView.hpp>

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
class Cylinder : public BasicGeometryEntity
{

  public:

    // Default constructor.
    Cylinder();

    // Constructor.
    Cylinder( const EntityId global_id, 
	      const int owner_rank, 
	      const int block_id,
	      const double length, 
	      const double radius,
	      const double centroid_x, 
	      const double centroid_y, 
	      const double centroid_z );

    // Destructor.
    ~Cylinder();

    //! Get the length of the cylinder.
    double length() const;

    //! Get the radius of the cylinder.
    double radius() const;

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

//! overload for printing cylinder
std::ostream& operator<< (std::ostream& os,const DataTransferKit::Cylinder& c); 

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_CYLINDER_HPP

//---------------------------------------------------------------------------//
// end DTK_Cylinder.hpp
//---------------------------------------------------------------------------//

