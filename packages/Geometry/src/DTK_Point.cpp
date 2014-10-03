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
 * \file DTK_Point.cpp
 * \author Stuart R. Slattery
 * \brief Point definition
 */
//---------------------------------------------------------------------------//

#include "DTK_Point.hpp"
#include "DTK_DBC.hpp"
#include "DTK_DataSerializer.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Default constructor.
Point::Point()
    : d_global_id( dtk_invalid_entity_id )
    , d_owner_rank( -1 )
    , d_coordinates( 0 )
    , d_byte_size( 0 )
{ /* ... */ }

//---------------------------------------------------------------------------//
// 1D Coordinate constructor.
Point::Point( const EntityId global_id, 
	      const int owner_rank,
	      const double x )
    : d_global_id( global_id )
    , d_owner_rank( owner_rank )
    , d_coordinates( 1 )
{
    d_coordinates[0] = x;
    d_byte_size = sizeof(EntityId) + sizeof(int) + sizeof(double);
}

//---------------------------------------------------------------------------//
// 2D Coordinate constructor.
Point::Point( const EntityId global_id, 
	      const int owner_rank, 
	      const double x, 
	      const double y )
    : d_global_id( global_id )
    , d_owner_rank( owner_rank )
    , d_coordinates( 2 )
{
    d_coordinates[0] = x;
    d_coordinates[1] = y;
    d_byte_size = sizeof(EntityId) + sizeof(int) + 2*sizeof(double);
}

//---------------------------------------------------------------------------//
// 3D Coordinate constructor.
Point::Point( const EntityId global_id, 
	      const int owner_rank,
	      const double x, 
	      const double y, 
	      const double z )
    : d_global_id( global_id )
    , d_owner_rank( owner_rank )
    , d_coordinates( 3 )
{
    d_coordinates[0] = x;
    d_coordinates[1] = y;
    d_coordinates[2] = z;
    d_byte_size = sizeof(EntityId) + sizeof(int) + 3*sizeof(double);
}

//---------------------------------------------------------------------------//
// Array constructor.
Point::Point( const EntityId global_id, 
	      const int owner_rank,
	      const Teuchos::Array<double>& coordinates )
    : d_global_id( global_id )
    , d_owner_rank( owner_rank )
    , d_coordinates( d_coordinates )
{
    d_byte_size = sizeof(EntityId) + sizeof(int) + 
		  (coordinates.size())*sizeof(double);
}

//---------------------------------------------------------------------------//
// Destructor.
Point::~Point()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the coordinates of the point.
void Point::getCoordinates( Teuchos::ArrayView<double>& coordinates ) const
{ 
    coordinates = d_coordinates(); 
}
//---------------------------------------------------------------------------//
// Return a string indicating the derived entity type.
std::string Point::entityType() const
{
    return std::string("DTK Point");
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
EntityId Point::id() const
{
    return d_global_id;
}
    
//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
int Point::ownerRank() const
{
    return d_owner_rank;
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
int Point::physicalDimension() const
{
    return d_coordinates.size();
}

//---------------------------------------------------------------------------//
// Return the parametric dimension of the entity.
int Point::parametricDimension() const
{
    return 0;
}

//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric
double Point::measure() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
void Point::centroid( const Teuchos::ArrayView<double>& centroid ) const
{
    getCoordinates( centroid );
}

//---------------------------------------------------------------------------//
// Return the axis-aligned bounding box around the entity.
void Point::boundingBox( Box& bounding_box ) const
{
    Teuchos::ArrayView<double> coordinates;
    getCoordinates( coordinates );
    box = Box( dtk_invalid_entity_type, d_owner_rank,
		coordinates[0], coordinates[1], coordinates[2], 
		coordinates[0], coordinates[1], coordinates[2] );
}

//---------------------------------------------------------------------------//
// Perform a safeguard check for mapping a point to the reference
void Point::safeguardMapToReferenceFrame(
    const Teuchos::ParameterList& parameters,
    const Teuchos::ArrayView<const double>& point,
    MappingStatus& status ) const
{
    bool not_implemented_for_point = true;
    DTK_INSIST( !not_implemented_for_point );
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the
void Point::mapToReferenceFrame( 
    const Teuchos::ParameterList& parameters,
    const Teuchos::ArrayView<const double>& point,
    const Teuchos::ArrayView<double>& reference_point,
    MappingStatus& status ) const
{
    bool not_implemented_for_point = true;
    DTK_INSIST( !not_implemented_for_point );
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of
bool Point::checkPointInclusion( 
    const Teuchos::ParameterList& parameters,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    bool not_implemented_for_point = true;
    DTK_INSIST( !not_implemented_for_point );
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
void Point::mapToPhysicalFrame( 
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& point ) const
{
    bool not_implemented_for_point = true;
    DTK_INSIST( !not_implemented_for_point );
}
     
//---------------------------------------------------------------------------//
// Get the size of the serialized entity in bytes.
std::size_t Point::byteSize() const
{
    DTK_REQUIRE( byte_size > 0 );
    return d_byte_size;
}

//---------------------------------------------------------------------------//
// Serialize the entity into a buffer.
void Point::serialize( const Teuchos::ArrayView<char>& buffer ) const
{
    DTK_REQUIRE( Teuchos::as<std::size_t>(buffer.size()) == d_byte_size );

    Teuchos::ArrayView<double> coordinates;
    getCoordinates( coordinates );
    std::size_t dimension = coordinates.size();

    DataSerializer serializer;
    serializer.setBuffer( buffer );
    serializer << d_global_id << d_owner_rank << dimension;

    Teuchos::ArrayView<double>::iterator coord_it;
    for ( coord_it = coordinates.begin(); 
	  coord_it != coordinates.end();
	  ++coord_it )
    {
	serializer << *coord_it;
    }
}

//---------------------------------------------------------------------------//
// Deserialize an entity from a buffer.
void Point::deserialize( const Teuchos::ArrayView<const char>& buffer )
{
    DTK_REQUIRE( Teuchos::as<std::size_t>(buffer.size()) == d_byte_size )

    Teuchos::ArrayView<char> buffer_nonconst(
	const_cast<char*>(buffer.getRawPtr()), buffer.size() );

    std::size_t dimension = 0;

    DataDeserializer deserializer;
    deserializer.setBuffer( buffer_nonconst );
    deserializer >> d_global_id >> d_owner_rank >> dimension;

    d_coordinates.resize( dimension );
    Teuchos::Array<double>::iterator coord_it;
    for ( coord_it = d_coordinates.begin(); 
	  coord_it != d_coordinates.end();
	  ++coord_it )
    {
	deserializer >> *coord_it;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Print the point description to an ostream.
 *
 * \return The ostream.
 */
std::ostream& operator<< (std::ostream& os,const DataTransferKit::Point& p)
{
    os << "Point: d_global_id=" << d_global_id 
       << ",d_owner_rank=" << d_owner_rank
       << ",d_coordinates=" << d_coordinates;

  return os;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_Point.cpp
//---------------------------------------------------------------------------//

