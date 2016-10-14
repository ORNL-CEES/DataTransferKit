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
 * \brief DTK_EntityLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief Forward and reverse local mappings for entities.
 */
//---------------------------------------------------------------------------//

#include <cmath>
#include <limits>

#include "DTK_DBC.hpp"
#include "DTK_EntityLocalMap.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
EntityLocalMap::EntityLocalMap() { /* ... */}

//---------------------------------------------------------------------------//
// Destructor.
EntityLocalMap::~EntityLocalMap() { /* ... */}

//---------------------------------------------------------------------------//
// Perform a safeguard check for mapping a point to the reference space
// of an entity using the given tolerance. Default implementation checks if
// the point is inside the bounding box of the entity.
bool EntityLocalMap::isSafeToMapToReferenceFrame(
    const Entity &entity, const Teuchos::ArrayView<const double> &point ) const
{
    // Get the bounding box of the entity.
    Teuchos::Tuple<double, 6> entity_box;
    entity.boundingBox( entity_box );

    // Check if the point is in the bounding box of the entity.
    double tolerance = 1.0e-6;
    int space_dim = entity.physicalDimension();
    bool in_x = true;
    if ( space_dim > 0 )
    {
        double x_tol = ( entity_box[3] - entity_box[0] ) * tolerance;
        in_x = ( ( point[0] >= ( entity_box[0] - x_tol ) ) &&
                 ( point[0] <= ( entity_box[3] + x_tol ) ) );
    }
    bool in_y = true;
    if ( space_dim > 1 )
    {
        double y_tol = ( entity_box[4] - entity_box[1] ) * tolerance;
        in_y = ( ( point[1] >= ( entity_box[1] - y_tol ) ) &&
                 ( point[1] <= ( entity_box[4] + y_tol ) ) );
    }
    bool in_z = true;
    if ( space_dim > 2 )
    {
        double z_tol = ( entity_box[5] - entity_box[2] ) * tolerance;
        in_z = ( ( point[2] >= ( entity_box[2] - z_tol ) ) &&
                 ( point[2] <= ( entity_box[5] + z_tol ) ) );
    }
    return ( in_x && in_y && in_z );
}

//---------------------------------------------------------------------------//
// Compute the normal on a face (3D) or edge (2D) at a given reference point.
void EntityLocalMap::normalAtReferencePoint(
    const Entity &entity, const Entity &parent_entity,
    const Teuchos::ArrayView<const double> &reference_point,
    const Teuchos::ArrayView<double> &normal ) const
{
    // Determine the reference dimension.
    int physical_dim = entity.physicalDimension();
    int ref_dim = physical_dim - 1;

    // Create a perturbation.
    double perturbation = std::sqrt( std::numeric_limits<double>::epsilon() );

    // 3D/face case.
    if ( 2 == ref_dim )
    {
        DTK_CHECK( 3 == reference_point.size() );
        DTK_CHECK( 3 == normal.size() );

        // Create extra points.
        Teuchos::Array<double> ref_p1( reference_point );
        Teuchos::Array<double> ref_p2( reference_point );

        // Apply a perturbation to the extra points.
        double p1_sign = 1.0;
        ref_p1[0] += perturbation;
        if ( !this->checkPointInclusion( entity, ref_p1() ) )
        {
            ref_p1[0] -= 2 * perturbation;
            p1_sign = -1.0;
        }
        double p2_sign = 1.0;
        ref_p2[1] += perturbation;
        if ( !this->checkPointInclusion( entity, ref_p2() ) )
        {
            ref_p2[1] -= 2 * perturbation;
            p2_sign = -1.0;
        }

        // Map the perturbed points to the physical frame.
        Teuchos::Array<double> p0( physical_dim );
        this->mapToPhysicalFrame( entity, reference_point(), p0() );
        Teuchos::Array<double> p1( physical_dim );
        this->mapToPhysicalFrame( entity, ref_p1(), p1() );
        Teuchos::Array<double> p2( physical_dim );
        this->mapToPhysicalFrame( entity, ref_p2(), p2() );

        // Compute the cross product of the tangents produced by the
        // perturbation.
        Teuchos::Array<double> tan1( physical_dim );
        Teuchos::Array<double> tan2( physical_dim );
        for ( int d = 0; d < physical_dim; ++d )
        {
            tan1[d] = p1_sign * ( p1[d] - p0[d] );
            tan2[d] = p2_sign * ( p2[d] - p0[d] );
        }
        normal[0] = tan1[1] * tan2[2] - tan1[2] * tan2[1];
        normal[1] = tan1[2] * tan2[0] - tan1[0] * tan2[2];
        normal[2] = tan1[0] * tan2[1] - tan1[1] * tan2[0];
    }

    // 2D/edge case.
    else if ( 1 == ref_dim )
    {
        DTK_CHECK( 2 == reference_point.size() );
        DTK_CHECK( 2 == normal.size() );

        // Create extra points.
        Teuchos::Array<double> ref_p1( reference_point );

        // Apply a perturbation to the extra points.
        double p1_sign = 1.0;
        ref_p1[0] += perturbation;
        if ( !this->checkPointInclusion( entity, ref_p1() ) )
        {
            ref_p1[0] -= 2 * perturbation;
            p1_sign = -1.0;
        }

        // Map the perturbed points to the physical frame.
        Teuchos::Array<double> p0( physical_dim );
        this->mapToPhysicalFrame( entity, reference_point(), p0() );
        Teuchos::Array<double> p1( physical_dim );
        this->mapToPhysicalFrame( entity, ref_p1(), p1() );

        // Compute the cross product of the tangents produced by the
        // perturbation.
        Teuchos::Array<double> tan( physical_dim );
        for ( int d = 0; d < physical_dim; ++d )
        {
            tan[d] = p1_sign * ( p1[d] - p0[d] );
        }
        normal[0] = -tan[0];
        normal[1] = tan[1];
    }

    // Normalize the normal vector.
    double norm = 0.0;
    for ( int d = 0; d < physical_dim; ++d )
    {
        norm += normal[d] * normal[d];
    }
    norm = std::sqrt( norm );
    for ( int d = 0; d < physical_dim; ++d )
    {
        normal[d] /= norm;
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_EntityLocalMap.cpp
//---------------------------------------------------------------------------//
