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
 * \brief DTK_LibmeshEntityLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief Forward and reverse local mappings for entities.
 */
//---------------------------------------------------------------------------//

#include <cassert>

#include "DTK_LibmeshEntityLocalMap.hpp"
#include "DTK_LibmeshHelpers.hpp"

#include <libmesh/fe_interface.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
LibmeshEntityLocalMap::LibmeshEntityLocalMap(
    const Teuchos::RCP<libMesh::MeshBase> &libmesh_mesh,
    const Teuchos::RCP<libMesh::System> &libmesh_system )
    : d_libmesh_mesh( libmesh_mesh )
    , d_libmesh_system( libmesh_system )
    , d_newton_tol( 1.0e-9 )
    , d_inclusion_tol( 1.0e-6 )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Set parameters for mapping.
void LibmeshEntityLocalMap::setParameters(
    const Teuchos::ParameterList &parameters )
{
    if ( parameters.isParameter( "Point Inclusion Tolerance" ) )
    {
        d_inclusion_tol = parameters.get<double>( "Point Inclusion Tolerance" );
    }
    if ( parameters.isParameter( "Newton Tolerance" ) )
    {
        d_newton_tol = parameters.get<double>( "Newton Tolerance" );
    }
}

//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric dimension (volume
// for a 3D entity, area for 2D, and length for 1D).
double
LibmeshEntityLocalMap::measure( const DataTransferKit::Entity &entity ) const
{
    if ( 0 == entity.topologicalDimension() )
    {
        return 0.0;
    }

    return LibmeshHelpers::extractGeom<libMesh::Elem>( entity )->volume();
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
void LibmeshEntityLocalMap::centroid(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<double> &centroid ) const
{
    libMesh::Point point;

    if ( 0 == entity.topologicalDimension() )
    {
        Teuchos::Ptr<libMesh::Point> node =
            LibmeshHelpers::extractGeom<libMesh::Node>( entity );
        point = *node;
    }
    else
    {
        point =
            LibmeshHelpers::extractGeom<libMesh::Elem>( entity )->centroid();
    }

    int space_dim = entity.physicalDimension();
    for ( int d = 0; d < space_dim; ++d )
    {
        centroid[d] = point( d );
    }
}

//---------------------------------------------------------------------------//
// Perform a safeguard check for mapping a point to the reference space
// of an entity using the given tolerance.
bool LibmeshEntityLocalMap::isSafeToMapToReferenceFrame(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &physical_point ) const
{
    int space_dim = entity.physicalDimension();
    int param_dim = 0;

    if ( 0 != entity.topologicalDimension() )
    {
        param_dim = LibmeshHelpers::extractGeom<libMesh::Elem>( entity )->dim();
    }
    else
    {
        return false;
    }

    if ( space_dim == param_dim )
    {
        // See if we are in the Cartesian bounding box.
        if ( EntityLocalMap::isSafeToMapToReferenceFrame( entity,
                                                          physical_point ) )
        {
            // If we are in the Cartesian bounding box see if we are 'close'
            // to the element according to libMesh.
            int space_dim = entity.physicalDimension();
            libMesh::Point lm_point;
            for ( int d = 0; d < space_dim; ++d )
            {
                lm_point( d ) = physical_point[d];
            }
            return LibmeshHelpers::extractGeom<libMesh::Elem>( entity )
                ->close_to_point( lm_point, d_inclusion_tol );
        }
        else
        {
            return false;
        }
    }
    else
    {
        // We currently do not natively support checks for mapping to
        // surfaces.
        bool not_implemented = true;
        assert( !not_implemented );
    }
    return false;
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the parameterized
// point.
bool LibmeshEntityLocalMap::mapToReferenceFrame(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &physical_point,
    const Teuchos::ArrayView<double> &reference_point ) const
{
    int space_dim = entity.physicalDimension();
    libMesh::Point lm_point;
    for ( int d = 0; d < space_dim; ++d )
    {
        lm_point( d ) = physical_point[d];
    }

    libMesh::Point lm_reference_point = libMesh::FEInterface::inverse_map(
        space_dim, d_libmesh_system->variable_type( 0 ),
        LibmeshHelpers::extractGeom<libMesh::Elem>( entity ).getRawPtr(),
        lm_point, d_newton_tol );

    for ( int d = 0; d < space_dim; ++d )
    {
        reference_point[d] = lm_reference_point( d );
    }

    return true;
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of an entity.
bool LibmeshEntityLocalMap::checkPointInclusion(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point ) const
{
    int space_dim = entity.physicalDimension();
    libMesh::Point lm_reference_point;
    for ( int d = 0; d < space_dim; ++d )
    {
        lm_reference_point( d ) = reference_point[d];
    }

    return libMesh::FEInterface::on_reference_element(
        lm_reference_point,
        LibmeshHelpers::extractGeom<libMesh::Elem>( entity )->type(),
        d_inclusion_tol );
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
void LibmeshEntityLocalMap::mapToPhysicalFrame(
    const DataTransferKit::Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    const Teuchos::ArrayView<double> &physical_point ) const
{
    int space_dim = entity.physicalDimension();
    libMesh::Point lm_reference_point;
    for ( int d = 0; d < space_dim; ++d )
    {
        lm_reference_point( d ) = reference_point[d];
    }

    libMesh::Point lm_point = libMesh::FEInterface::map(
        space_dim, d_libmesh_system->variable_type( 0 ),
        LibmeshHelpers::extractGeom<libMesh::Elem>( entity ).getRawPtr(),
        lm_reference_point );

    for ( int d = 0; d < space_dim; ++d )
    {
        physical_point[d] = lm_point( d );
    }
}

//---------------------------------------------------------------------------//
// Compute the normal on a face (3D) or edge (2D) at a given reference point.
void LibmeshEntityLocalMap::normalAtReferencePoint(
    const DataTransferKit::Entity &entity,
    const DataTransferKit::Entity &parent_entity,
    const Teuchos::ArrayView<const double> &reference_point,
    const Teuchos::ArrayView<double> &normal ) const
{
    DataTransferKit::EntityLocalMap::normalAtReferencePoint(
        entity, parent_entity, reference_point, normal );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_LibmeshEntityLocalMap.cpp
//---------------------------------------------------------------------------//
