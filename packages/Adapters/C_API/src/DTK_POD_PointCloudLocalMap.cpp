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
 * \brief DTK_POD_PointCloudLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief Forward and reverse local mappings for entities.
 */
//---------------------------------------------------------------------------//

#include "DTK_POD_PointCloudLocalMap.hpp"
#include "DTK_POD_PointCloudEntity.hpp"
#include "DTK_POD_PointCloudEntityImpl.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Set parameters for mapping.
void POD_PointCloudLocalMap::setParameters(
    const Teuchos::ParameterList& parameters )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric dimension (volume
// for a 3D entity, area for 2D, and length for 1D).
double POD_PointCloudLocalMap::measure( const Entity& entity ) const
{
    return 0.0;
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
void POD_PointCloudLocalMap::centroid(
    const Entity& entity, const Teuchos::ArrayView<double>& centroid ) const
{
    for ( int d = 0; d < entity.physicalDimension(); ++d )
    {
        centroid[d] =
            Teuchos::rcp_dynamic_cast<POD_PointCloudEntityImpl>(
                entity.extraData())->coord(d);
    }
}

//---------------------------------------------------------------------------//
// Perform a safeguard check for mapping a point to the reference space
// of an entity using the given tolerance.
bool POD_PointCloudLocalMap::isSafeToMapToReferenceFrame(
    const Entity& entity,
    const Teuchos::ArrayView<const double>& physical_point ) const
{
    // Points have no reference frame.
    return false;
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the parameterized
// point.
bool POD_PointCloudLocalMap::mapToReferenceFrame(
    const Entity& entity,
    const Teuchos::ArrayView<const double>& physical_point,
    const Teuchos::ArrayView<double>& reference_point ) const
{
    // Points have no reference frame.
    return false;
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of an entity.
bool POD_PointCloudLocalMap::checkPointInclusion(
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    // Always false.
    return false;
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
void POD_PointCloudLocalMap::mapToPhysicalFrame(
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& physical_point ) const
{
    // Do nothing.
}

//---------------------------------------------------------------------------//
// Compute the normal on a face (3D) or edge (2D) at a given reference point.
void POD_PointCloudLocalMap::normalAtReferencePoint(
    const Entity& entity,
    const Entity& parent_entity,
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& normal ) const
{
    // Do nothing.
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_POD_PointCloudLocalMap.cpp
//---------------------------------------------------------------------------//
