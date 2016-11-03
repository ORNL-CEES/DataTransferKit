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
 * \brief DTK_Jacobian.cpp
 * \author Stuart R. Slattery
 * \brief Geometric entity set interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_Jacobian.hpp"
#include "DTK_BasicEntityPredicates.hpp"
#include "DTK_DBC.hpp"

#include <Teuchos_CommHelpers.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
Jacobian::Jacobian( const Teuchos::RCP<FunctionSpace> &function_space )
    : d_function_space( function_space )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Destructor.
Jacobian::~Jacobian() { /* ... */}

//---------------------------------------------------------------------------//
// Get the Jacobian for a reference point in an entity
void Jacobian::jacobian(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    Teuchos::Array<Teuchos::Array<double>> &jac ) const
{
    auto local_map = d_function_space->localMap();
    auto set = d_function_space->entitySet();
    auto shape_function = d_function_space->shapeFunction();
    int space_dim = entity.topologicalDimension();

    // Get the support ids of the entity.
    Teuchos::Array<SupportId> support_ids;
    shape_function->entitySupportIds( entity, support_ids );
    int cardinality = support_ids.size();

    // Get physical coordinates of support
    Teuchos::Array<Teuchos::Array<double>> support_coordinates( cardinality );
    for ( int ni = 0; ni < cardinality; ni++ )
    {
        DataTransferKit::Entity support;
        set->getEntity( support_ids[ni], 0, support );

        support_coordinates[ni].resize( space_dim );
        local_map->centroid( support, support_coordinates[ni]() );
    }

    // Evaluate the gradients of shape functions.
    Teuchos::Array<Teuchos::Array<double>> shape_grads;
    shape_function->evaluateGradient( entity, reference_point, shape_grads );

    // Compute Jacobian
    jac.resize( space_dim );
    for ( int i = 0; i < space_dim; i++ )
    {
        jac[i].resize( space_dim );
        for ( int j = 0; j < space_dim; j++ )
        {
            jac[i][j] = 0.0;
            for ( int ni = 0; ni < cardinality; ++ni )
                jac[i][j] += support_coordinates[ni][i] * shape_grads[ni][j];
        }
    }
}

// Get the determinant of the Jacobian for a reference point in an entity
double Jacobian::jacobian_determinant(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point ) const
{
    int space_dim = entity.topologicalDimension();
    DTK_REQUIRE( space_dim == 2 || space_dim == 3 );

    Teuchos::Array<Teuchos::Array<double>> J;
    jacobian( entity, reference_point, J );

    // Compute derminant.
    double det = 0.0;
    if ( space_dim == 2 )
    {
        det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    }
    else if ( space_dim == 3 )
    {
        det = J[0][0] * J[1][1] * J[2][2] + J[2][0] * J[0][1] * J[1][2] +
              J[0][2] * J[1][0] * J[2][1] - J[0][2] * J[1][1] * J[2][0] -
              J[0][0] * J[1][2] * J[2][1] - J[2][2] * J[1][0] * J[0][1];
    }

    return det;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_Jacobian.cpp
//---------------------------------------------------------------------------//
