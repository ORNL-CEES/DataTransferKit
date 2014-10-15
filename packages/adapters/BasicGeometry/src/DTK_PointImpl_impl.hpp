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
 * \file DTK_Point_impl.hpp
 * \author Stuart R. Slattery
 * \brief Point definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_POINTIMPL_IMPL_HPP
#define DTK_POINTIMPL_IMPL_HPP

#include <limits>

#include "DTK_DBC.hpp"
#include "DTK_DataSerializer.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Default constructor.
template<int DIM>
PointImpl<DIM>::PointImpl()
    : d_global_id( dtk_invalid_entity_id )
    , d_owner_rank( -1 )
    , d_coordinates( 0 )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Array constructor.
template<int DIM>
PointImpl<DIM>::PointImpl( const EntityId global_id, 
			   const int owner_rank,
			   const Teuchos::Array<double>& coordinates )
    : d_global_id( global_id )
    , d_owner_rank( owner_rank )
    , d_coordinates( coordinates )
{
    DTK_REQUIRE( DIM == Teuchos::as<int>(d_coordinates.size()) );
}

//---------------------------------------------------------------------------//
// Destructor.
template<int DIM>
PointImpl<DIM>::~PointImpl()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the coordinates of the point.
template<int DIM>
void PointImpl<DIM>::getCoordinates( 
    Teuchos::ArrayView<const double>& coordinates ) const
{ 
    coordinates = d_coordinates(); 
}
//---------------------------------------------------------------------------//
// Return a string indicating the derived entity type.
template<int DIM>
std::string PointImpl<DIM>::name() const
{
    return std::string("DTK Point");
}

//---------------------------------------------------------------------------//
// Get the entity type.
template<int DIM>
EntityType PointImpl<DIM>::entityType() const
{
    return NODE;
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
template<int DIM>
EntityId PointImpl<DIM>::id() const
{
    return d_global_id;
}
    
//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
template<int DIM>
int PointImpl<DIM>::ownerRank() const
{
    return d_owner_rank;
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
template<int DIM>
int PointImpl<DIM>::physicalDimension() const
{
    return DIM;
}

//---------------------------------------------------------------------------//
// Return the parametric dimension of the entity.
template<int DIM>
int PointImpl<DIM>::parametricDimension() const
{
    return 0;
}

//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric
template<int DIM>
double PointImpl<DIM>::measure() const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
template<int DIM>
void PointImpl<DIM>::centroid( Teuchos::ArrayView<const double>& centroid ) const
{
    getCoordinates( centroid );
}

//---------------------------------------------------------------------------//
// Return the axis-aligned bounding box around the entity. 1D specialization.
template<>
void PointImpl<1>::boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    Teuchos::ArrayView<const double> coordinates;
    getCoordinates( coordinates );
    double max = std::numeric_limits<double>::max();
    bounds = Teuchos::tuple( coordinates[0], -max, -max,
			     coordinates[0], max, max );
}

//---------------------------------------------------------------------------//
// Return the axis-aligned bounding box around the entity. 2D specialization.
template<>
void PointImpl<2>::boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    Teuchos::ArrayView<const double> coordinates;
    getCoordinates( coordinates );
    double max = std::numeric_limits<double>::max();
    bounds = Teuchos::tuple( coordinates[0], coordinates[1], -max,
			     coordinates[0], coordinates[1], max );
}

//---------------------------------------------------------------------------//
// Return the axis-aligned bounding box around the entity. 3D specialization.
template<>
void PointImpl<3>::boundingBox( Teuchos::Tuple<double,6>& bounds ) const
{
    Teuchos::ArrayView<const double> coordinates;
    getCoordinates( coordinates );
    bounds = Teuchos::tuple( coordinates[0], coordinates[1], coordinates[2], 
			     coordinates[0], coordinates[1], coordinates[2] );
}

//---------------------------------------------------------------------------//
// Perform a safeguard check for mapping a point to the reference
template<int DIM>
void PointImpl<DIM>::safeguardMapToReferenceFrame(
    const Teuchos::ParameterList& parameters,
    const Teuchos::ArrayView<const double>& point,
    MappingStatus& status ) const
{
    bool not_implemented_for_point = true;
    DTK_INSIST( !not_implemented_for_point );
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the
template<int DIM>
void PointImpl<DIM>::mapToReferenceFrame( 
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
template<int DIM>
bool PointImpl<DIM>::checkPointInclusion( 
    const Teuchos::ParameterList& parameters,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    bool not_implemented_for_point = true;
    DTK_INSIST( !not_implemented_for_point );
    return false;
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
template<int DIM>
void PointImpl<DIM>::mapToPhysicalFrame( 
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& point ) const
{
    bool not_implemented_for_point = true;
    DTK_INSIST( !not_implemented_for_point );
}
     
//---------------------------------------------------------------------------//
// Serialize the entity into a buffer.
template<int DIM>
void PointImpl<DIM>::serialize( const Teuchos::ArrayView<char>& buffer ) const
{
    DTK_REQUIRE( Teuchos::as<std::size_t>(buffer.size()) >= 
		 PointImpl<DIM>::byteSize() );

    Teuchos::ArrayView<const double> coordinates;
    getCoordinates( coordinates );

    DataSerializer serializer;
    serializer.setBuffer( buffer );
    serializer << d_global_id << d_owner_rank;

    Teuchos::ArrayView<const double>::const_iterator coord_it;
    for ( coord_it = coordinates.begin(); 
	  coord_it != coordinates.end();
	  ++coord_it )
    {
	serializer << *coord_it;
    }
}

//---------------------------------------------------------------------------//
// Deserialize an entity from a buffer.
template<int DIM>
void PointImpl<DIM>::deserialize( const Teuchos::ArrayView<const char>& buffer )
{
    DTK_REQUIRE( Teuchos::as<std::size_t>(buffer.size()) >=
		 sizeof(EntityId) + sizeof(int) + DIM*sizeof(double) );
    Teuchos::ArrayView<char> buffer_nonconst(
	const_cast<char*>(buffer.getRawPtr()), buffer.size() );

    DataDeserializer deserializer;
    deserializer.setBuffer( buffer_nonconst );
    deserializer >> d_global_id >> d_owner_rank;

    d_coordinates.resize( DIM );
    Teuchos::Array<double>::iterator coord_it;
    for ( coord_it = d_coordinates.begin(); 
	  coord_it != d_coordinates.end();
	  ++coord_it )
    {
	deserializer >> *coord_it;
    }
}

//---------------------------------------------------------------------------//
// Static Members.
//---------------------------------------------------------------------------//
// Byte size of the point.
template<int DIM>
std::size_t 
PointImpl<DIM>::d_byte_size = sizeof(EntityId) + sizeof(int) + DIM*sizeof(double);

//---------------------------------------------------------------------------//
// Get the byte size of the point.
template<int DIM>
std::size_t PointImpl<DIM>::byteSize()
{
    return d_byte_size;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_POINTIMPL_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_PointImpl_impl.hpp
//---------------------------------------------------------------------------//

