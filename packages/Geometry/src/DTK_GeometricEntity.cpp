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
 * \brief DTK_GeometricEntity.cpp
 * \author Stuart R. Slattery
 * \brief Geometric entity interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_GeometricEntity.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
GeometricEntity::GeometricEntity()
{ /* ... */ }

//---------------------------------------------------------------------------//
//brief Destructor.
GeometricEntity::~GeometricEntity()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
EntityId GeometricEntity::id() const
{ 
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return -1;
}
    
//---------------------------------------------------------------------------//
// Get the native parallel rank of the entity.
int GeometricEntity::parallelRank() const
{ 
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return -1;
}
//---------------------------------------------------------------------------//
// Return the spatial dimension of the entity.
int GeometricEntity::dimension() const
{ 
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return -1;
}

//---------------------------------------------------------------------------//
// Return the entity measure (volume for a 3D entity, area for 2D, and length
// for 1D). 
double GeometricEntity::measure() const
{ 
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return -1.0;
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
Teuchos::Array<double> centroid() const
{ 
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return Teuchos::Array<double>(0);
}

//---------------------------------------------------------------------------//
// Return the axis-aligned bounding box around the entity.
Box GeometricEntity::boundingBox() const
{ 
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return BoundingBox();
}

//---------------------------------------------------------------------------//
// Perform a safeguard check for mapping a point to the reference space of an
// entity using the given tolerance. 
bool GeometricEntity::isSafeToMapToReferenceFrame(
    const Teuchos::ArrayView<double>& point,
    const double tolerance ) const
{ 
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return false;
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the
// parameterized point. 
void GeometricEntity::mapToReferenceFrame( 
    const Teuchos::ArrayView<double>& point,
    const double tolerance,
    Teuchos::Array<double>& reference_point ) const
{ 
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of an entity.
bool GeometricEntity::checkPointInclusion( 
    const Teuchos::Array<double>& reference_point ) const
{ 
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
    return false;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_GeometricEntity.cpp
//---------------------------------------------------------------------------//
