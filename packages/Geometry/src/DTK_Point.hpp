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

#include "DTK_GeometryTypes.hpp"
#include "DTK_GeometricEntity.hpp"

#include <Teuchos_ArrayView.hpp>

#include <iostream>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class Point
 * \brief Point container.
 *
 * Users can subclass this function for more efficient access to and storage
 * of coordinates. If a set of point coordinates already existed, a subclass
 * could be used to point to those coordinates through this interface. This
 * interface does assume access to contiguous storage of (x,y,z) coordinates
 * for a given point. 
 */
//---------------------------------------------------------------------------//
class Point : public GeometricEntity
{
  public:

    // Default constructor.
    Point();

    // 1D Coordinate constructor.
    Point( const EntityId global_id, 
	   const int owner_rank,
	   const double x );

    // 2D Coordinate constructor.
    Point( const EntityId global_id, 
	   const int owner_rank, 
	   const double x, 
	   const double y );

    // 3D Coordinate constructor.
    Point( const EntityId global_id, 
	   const int owner_rank,
	   const double x, 
	   const double y, 
	   const double z );

    // Array constructor.
    Point( const EntityId global_id, 
	   const int owner_rank,
	   const Teuchos::Array<double>& coordinates );

    // Destructor.
    ~Point();

    //@{
    //! Coordinate access functions.
    // Get the coordinates of the point.
    virtual void getCoordinates( Teuchos::ArrayView<double>& coordinates ) const;
    //@}

    //@{ 
    //! GeometricEntity implementation.
    // Return a string indicating the derived entity type.
    std::string entityType() const;

    // Get the unique global identifier for the entity.
    EntityId id() const;
    
    // Get the parallel rank that owns the entity.
    int ownerRank() const;

    // Return the physical dimension of the entity.
    virtual int physicalDimension() const;

    // Return the parametric dimension of the entity.
    int parametricDimension() const;

    // Return the entity measure with respect to the parameteric
    double measure() const;

    // Return the centroid of the entity.
    void centroid( const Teuchos::ArrayView<double>& centroid ) const;

    // Return the axis-aligned bounding box around the entity.
    void boundingPoint( Point& bounding_box ) const;

    // Perform a safeguard check for mapping a point to the reference
    void safeguardMapToReferenceFrame(
	const Teuchos::ParameterList& parameters,
	const Teuchos::ArrayView<const double>& point,
	MappingStatus& status ) const;

    // Map a point to the reference space of an entity. Return the
    void mapToReferenceFrame( 
	const Teuchos::ParameterList& parameters,
	const Teuchos::ArrayView<const double>& point,
	const Teuchos::ArrayView<double>& reference_point,
	MappingStatus& status ) const;

    // Determine if a reference point is in the parameterized space of
    bool checkPointInclusion( 
	const Teuchos::ParameterList& parameters,
	const Teuchos::ArrayView<const double>& reference_point ) const;

    // Map a reference point to the physical space of an entity.
    void mapToPhysicalFrame( 
	const Teuchos::ArrayView<const double>& reference_point,
	const Teuchos::ArrayView<double>& point ) const;
     
    // Get the size of the serialized entity in bytes.
    std::size_t byteSize() const;

    // Serialize the entity into a buffer.
    void serialize( const Teuchos::ArrayView<char>& buffer ) const;

    // Deserialize an entity from a buffer.
    void deserialize( const Teuchos::ArrayView<const char>& buffer );
    //@}

  protected:

    // Global id.
    EntityId d_global_id;

    // Owning parallel rank.
    int d_owner_rank;

  private:

    // Coordinates.
    Teuchos::Array<double> d_coordinates;

    // Packed size in bytes.
    std::size_t d_byte_size;
};

//---------------------------------------------------------------------------//
//! overload for printing Point
std::ostream& operator<< (std::ostream& os,const DataTransferKit::Point& p); 

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_POINT_HPP

//---------------------------------------------------------------------------//
// end DTK_Point.hpp
//---------------------------------------------------------------------------//

