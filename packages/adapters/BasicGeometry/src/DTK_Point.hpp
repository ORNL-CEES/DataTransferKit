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
 * \file DTK_Point.hpp
 * \author Stuart R. Slattery
 * \brief Point declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_POINT_HPP
#define DTK_POINT_HPP

#include <iostream>

#include "DTK_PointImpl.hpp"
#include "DTK_Types.hpp"
#include "DTK_Entity.hpp"
#include "DTK_Box.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class Point
 * \brief Point container declaration.
 */
//---------------------------------------------------------------------------//
class Point : public Entity
{
  public:

    // Default constructor.
    Point();

    // Array constructor.
    Point( const EntityId global_id, 
	   const int owner_rank,
	   const Teuchos::Array<double>& coordinates,
	   const bool on_surface = false,
	   const Teuchos::Array<int>& block_ids = Teuchos::Array<int>(0),
	   const Teuchos::Array<int>& boundary_ids = Teuchos::Array<int>(0) );

    // Destructor.
    ~Point();

    //@{
    //! Coordinate access functions.
    // Get the coordinates of the point.
    void getCoordinates( const Teuchos::ArrayView<double>& coordinates ) const;
    //@}

    // Return the entity measure with respect to the parameteric
    double measure() const;

    // Compute the centroid of the entity.
    void centroid( const Teuchos::ArrayView<double>& centroid ) const;

    // (Safeguard the reverse map) Perform a safeguard check for mapping a
    // point to the reference space of an entity using the given tolerance.
    bool isSafeToMapToReferenceFrame(
	const Teuchos::ArrayView<const double>& point ) const;

    // (Reverse Map) Map a point to the reference space of an entity. Return
    // the parameterized point.
    bool mapToReferenceFrame( 
	const Teuchos::ArrayView<const double>& point,
	const Teuchos::ArrayView<double>& reference_point ) const;

    // Determine if a reference point is in the parameterized space of an
    // entity.
    bool checkPointInclusion( 
	const double tolerance,
	const Teuchos::ArrayView<const double>& reference_point ) const;

    // (Forward Map) Map a reference point to the physical space of an entity.
    void mapToPhysicalFrame( 
	const Teuchos::ArrayView<const double>& reference_point,
	const Teuchos::ArrayView<double>& point ) const;

  private:

    // Point implementation.
    Teuchos::RCP<PointImpl > d_point_impl;
};

//---------------------------------------------------------------------------//
//! overload for printing Point
template<int DIM>
std::ostream& operator<< (std::ostream& os,const DataTransferKit::Point& p); 

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_POINT_HPP

//---------------------------------------------------------------------------//
// end DTK_Point.hpp
//---------------------------------------------------------------------------//

