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

#include "DTK_PointImpl.hpp"
#include "DTK_Types.hpp"
#include "DTK_Entity.hpp"
#include "DTK_Box.hpp"
#include "DTK_AbstractObjectRegistry.hpp"

#include <Teuchos_ArrayView.hpp>

#include <iostream>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class Point
 * \brief Point container declaration.
 */
//---------------------------------------------------------------------------//
template<int DIM>
class Point : public Entity
{
  public:

    // Default constructor.
    Point();

    // Array constructor.
    Point( const EntityId global_id, 
	   const int owner_rank,
	   const Teuchos::Array<double>& coordinates );

    // Destructor.
    ~Point();

    //@{
    //! Coordinate access functions.
    // Get the coordinates of the point.
    void getCoordinates( Teuchos::ArrayView<const double>& coordinates ) const;
    //@}
    //@}

    // Get the byte size for the point.
    static std::size_t byteSize();

  private:

    // Point implementation.
    Teuchos::RCP<PointImpl<DIM> > d_point_impl;
};

//---------------------------------------------------------------------------//
//! overload for printing Point
template<int DIM>
std::ostream& operator<< (std::ostream& os,const DataTransferKit::Point<DIM>& p); 

//---------------------------------------------------------------------------//
// AbstractObjectRegistrationPolicy implementation.
//---------------------------------------------------------------------------//
template<int DIM>
class AbstractObjectRegistrationPolicy<Point<DIM> >
{
  public:

    //! Base class type.
    typedef Point<DIM> object_type;

    /*!
     * \brief Register a derived class with a base class.
     */
    static void registerDerivedClassWithBaseClass()
    {
	// Register the constructor with the base class
	// AbstractBuildableObject interface.
	Entity::setDerivedClassFactory<Point<DIM> >();

	// Register the byte size with the base class
	// AbstractSerializableObject interface.
	Entity::setDerivedClassByteSize( Point<DIM>::byteSize() );
    }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_Point_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_POINT_HPP

//---------------------------------------------------------------------------//
// end DTK_Point.hpp
//---------------------------------------------------------------------------//

