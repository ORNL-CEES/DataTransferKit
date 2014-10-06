//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \brief Box declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BOX_HPP
#define DTK_BOX_HPP

#include "DTK_GeometricEntity.hpp"
#include "DTK_DerivedObjectRegistry.hpp"
#include "DTK_BoxImpl.hpp"

#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_SerializationTraits.hpp>

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
class Box : public GeometricEntity
{
  public:

    // Default constructor.
    Box();

    // Constructor.
    Box( const EntityId global_id, const int owner_rank,
	 const double x_min, const double y_min, const double z_min,
	 const double x_max, const double y_max, const double z_max );

    // Tuple constructor.
    Box( const EntityId global_id,
	 const int owner_rank, 
	 const Teuchos::Tuple<double,6>& bounds );

    // Destructor.
    ~Box();

    // Get the boundaries of the box.
    Teuchos::Tuple<double,6> getBounds() const
    { return d_box_impl->getBounds(); }

    // Static function to check for box intersection but not perform it.
    static bool checkForIntersection( const Box& box_A,
				      const Box& box_B );

    // Static function for box intersection.
    static bool intersectBoxes( const Box& box_A,
				const Box& box_B,
				Box& box_intersection );

    // Static function for box union
    static void uniteBoxes( const Box& box_A,
			    const Box& box_B,
			    Box& box_union );

    // Compound assignment overload.
    Box& operator+=(const Box& rhs);

    // Get the byte size for the box.
    static std::size_t byteSize();

  private:

    // Pointer to the implementation.
    Teuchos::RCP<BoxImpl> d_box_impl;
};

//---------------------------------------------------------------------------//
//! Addition operator overload. Adding two boxes together is the same as
//! computing their union.
Box operator+( const Box& box_1, const Box& box_2 );

//---------------------------------------------------------------------------//
//! overload for printing box
std::ostream& operator<< (std::ostream& os,const DataTransferKit::Box& b); 

//---------------------------------------------------------------------------//
// DerivedObjectRegistrationPolicy implementation.
//---------------------------------------------------------------------------//
template<>
class DerivedObjectRegistrationPolicy<Box>
{
  public:

    //! Base class type.
    typedef Box object_type;

    /*!
     * \brief Register a derived class with a base class.
     */
    static void registerDerivedClassWithBaseClass()
    {
	// Register the constructor with the base class
	// AbstractBuildableObject interface.
	GeometricEntity::setDerivedClassFactory<Box>();

	// Register the byte size with the base class
	// AbstractSerializableObject interface.
	GeometricEntity::setDerivedClassByteSize( Box::byteSize() );
    }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Teuchos::SerializationTraits implementation.
//---------------------------------------------------------------------------//

namespace Teuchos
{
template<typename Ordinal>
class SerializationTraits<Ordinal,DataTransferKit::Box> 
{
  public:

    static const bool supportsDirectSerialization = false;

    static Ordinal fromCountToIndirectBytes( const Ordinal count, 
					     const DataTransferKit::Box buffer[] ) 
    { 
	return count * DataTransferKit::Box::byteSize();
    }

    static void serialize( const Ordinal count, 
			   const DataTransferKit::Box buffer[], 
			   const Ordinal bytes, 
			   char charBuffer[] )
    {
	DTK_REQUIRE( fromCountToIndirectBytes(count,buffer) == bytes );
	std::size_t box_size = DataTransferKit::Box::byteSize();
	char* buffer_pos = &charBuffer[0];
	for ( int n = 0; n < count; ++n )
	{
	    Teuchos::ArrayView<char> buffer_view( buffer_pos, box_size );
	    buffer[n].serialize( buffer_view );
	    buffer_pos += box_size;
	}
	DTK_CHECK( &charBuffer[0] + bytes == buffer_pos);
    }

    static Ordinal fromIndirectBytesToCount( const Ordinal bytes, 
					     const char charBuffer[] ) 
    { 
	return bytes / DataTransferKit::Box::byteSize();
    }

    static void deserialize( const Ordinal bytes, 
			     const char charBuffer[], 
			     const Ordinal count, 
			     DataTransferKit::Box buffer[] )
    { 
	DTK_REQUIRE( fromIndirectBytesToCount(bytes,charBuffer) == count );
	std::size_t box_size = DataTransferKit::Box::byteSize();
	char* buffer_pos = const_cast<char*>(&charBuffer[0]);
	for ( int n = 0; n < count; ++n )
	{
	    Teuchos::ArrayView<char> buffer_view( buffer_pos, box_size );
	    buffer[n].deserialize( buffer_view );
	    buffer_pos += box_size;
	}
	DTK_CHECK( &charBuffer[0] + bytes == buffer_pos );
    }
};
} // end namespace Teuchos

//---------------------------------------------------------------------------//

#endif // end DTK_BOX_HPP

//---------------------------------------------------------------------------//
// end DTK_Box.hpp
//---------------------------------------------------------------------------//

