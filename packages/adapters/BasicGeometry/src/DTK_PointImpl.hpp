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

#ifndef DTK_POINTIMPL_HPP
#define DTK_POINTIMPL_HPP

#include "DTK_Types.hpp"
#include "DTK_EntityImpl.hpp"

#include <Teuchos_ArrayView.hpp>

#include <iostream>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class Point
 * \brief Point container implementation.
 *
 * Users can subclass this class for more efficient access to and storage of
 * point coordinates while getting the rest of the Entity
 * implementation for free. For example, if a set of point coordinates already
 * existed, a subclass could be used to point to those coordinates through
 * this interface. This interface does assume access to contiguous storage of
 * (x,y,z) coordinates for a given point such that a pointer to these coordinates
 * can be accessed with a Teuchos::ArrayView. Serializing a point accesses
 * coordinates through this interface. Deserializing a point constructs
 * this class directly instead of potential subclasses. If in the future we
 * want to directly deserialize to a subclass of Point through the geometric
 * entity interface we can make the name and serialization functions
 * virtual to permit the subclass to override them. A user could also just
 * directly subclass Entity for their particular point type.
 */
//---------------------------------------------------------------------------//
class PointImpl : public EntityImpl
{
  public:

    // Default constructor.
    PointImpl();

    // Array constructor.
    PointImpl( const EntityId global_id, 
	       const int owner_rank,
	       const Teuchos::Array<double>& coordinates );

    // Destructor.
    ~PointImpl();

    //@{
    //! Coordinate access functions.
    // Get the coordinates of the point.
    virtual void 
    getCoordinates( Teuchos::ArrayView<const double>& coordinates ) const;
    //@}

    // Get the entity type.
    EntityType entityType() const;

    // Get the unique global identifier for the entity.
    EntityId id() const;
    
    // Get the parallel rank that owns the entity.
    int ownerRank() const;

    // Return the physical dimension of the entity.
    virtual int physicalDimension() const;
     
  protected:

    // Global id.
    EntityId d_global_id;

    // Owning parallel rank.
    int d_owner_rank;

  private:

    // Coordinates.
    Teuchos::Array<double> d_coordinates;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_POINTIMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_PointImpl.hpp
//---------------------------------------------------------------------------//

