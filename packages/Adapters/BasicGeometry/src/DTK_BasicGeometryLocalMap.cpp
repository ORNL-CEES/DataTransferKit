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
 * \brief DTK_EntityLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief Forward and reverse local mappings for entities.
 */
//---------------------------------------------------------------------------//

#include "DTK_BasicGeometryLocalMap.hpp"
#include "DTK_BasicGeometryEntity.hpp"
#include "DTK_BasicGeometryExtraData.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
BasicGeometryLocalMap::BasicGeometryLocalMap()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric dimension (volume
// for a 3D entity, area for 2D, and length for 1D). 
double BasicGeometryLocalMap::measure( const Entity& entity ) const
{
    return Teuchos::rcp_dynamic_cast<BasicGeometryExtraData>(entity.extraData()
	)->implementationConstPtr()->measure();
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
void BasicGeometryLocalMap::centroid( 
    const Entity& entity, const Teuchos::ArrayView<double>& centroid ) const
{ 
    Teuchos::rcp_dynamic_cast<BasicGeometryExtraData>(entity.extraData()
	)->implementationConstPtr()->centroid(centroid);
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the parameterized point.
bool BasicGeometryLocalMap::mapToReferenceFrame( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& point,
    const Teuchos::ArrayView<double>& reference_point,
    const Teuchos::RCP<MappingStatus>& status ) const
{
    return Teuchos::rcp_dynamic_cast<BasicGeometryExtraData>(entity.extraData()
	)->implementationConstPtr()->mapToReferenceFrame(point,reference_point);
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of an entity.
bool BasicGeometryLocalMap::checkPointInclusion( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    double tolerance = 1.0e-6;
    if ( Teuchos::nonnull(this->b_parameters) )
    {
	if ( this->b_parameters->isParameter("Point Inclusion Tolerance") )
	{
	    tolerance = 
		this->b_parameters->get<double>("Point Inclusion Tolerance");
	}
    }
    return Teuchos::rcp_dynamic_cast<BasicGeometryExtraData>(entity.extraData()
	)->implementationConstPtr()->checkPointInclusion(tolerance,reference_point);
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
void BasicGeometryLocalMap::mapToPhysicalFrame( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& point ) const
{
    Teuchos::rcp_dynamic_cast<BasicGeometryExtraData>(entity.extraData()
	)->implementationConstPtr()->mapToPhysicalFrame(reference_point,point);
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_BasicGeometryLocalMap.cpp
//---------------------------------------------------------------------------//
