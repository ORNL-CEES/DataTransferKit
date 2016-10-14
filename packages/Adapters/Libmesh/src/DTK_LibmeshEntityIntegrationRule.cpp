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
 * \brief DTK_LibmeshEntityIntegrationRule.cpp
 * \author Stuart R. Slattery
 * \brief libMesh integration rule implementation.
 */
//---------------------------------------------------------------------------//

#include "DTK_LibmeshEntityIntegrationRule.hpp"
#include "DTK_LibmeshHelpers.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
LibmeshEntityIntegrationRule::LibmeshEntityIntegrationRule(
    const libMesh::QuadratureType quadrature_type )
    : d_quad_type( quadrature_type )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Given an entity and an integration order, get its integration rule.
void LibmeshEntityIntegrationRule::getIntegrationRule(
    const Entity &entity, const int order,
    Teuchos::Array<Teuchos::Array<double>> &reference_points,
    Teuchos::Array<double> &weights ) const
{
    // Create a libmesh quadrature rule. Use Gauss quadrature as the default.
    libMesh::UniquePtr<libMesh::QBase> libmesh_quadrature =
        libMesh::QBase::build( d_quad_type, entity.topologicalDimension(),
                               static_cast<libMesh::Order>( order ) );

    // Initialize the quadrature rule for the entity.
    Teuchos::Ptr<libMesh::Elem> libmesh_elem =
        LibmeshHelpers::extractGeom<libMesh::Elem>( entity );
    libmesh_quadrature->init( libmesh_elem->type() );

    // Extract the data for the quadrature rule.
    int num_points = libmesh_quadrature->n_points();
    int quad_dim = libmesh_quadrature->get_dim();
    reference_points.resize( num_points );
    weights.resize( num_points );
    for ( int p = 0; p < num_points; ++p )
    {
        weights[p] = libmesh_quadrature->w( p );
        libMesh::Point qp = libmesh_quadrature->qp( p );
        reference_points[p].resize( quad_dim );
        for ( int d = 0; d < quad_dim; ++d )
        {
            reference_points[p][d] = qp( d );
        }
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_LibmeshEntityIntegrationRule.hpp
//---------------------------------------------------------------------------//
