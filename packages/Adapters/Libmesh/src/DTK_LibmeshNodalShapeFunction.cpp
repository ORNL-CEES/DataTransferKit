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
 * \brief DTK_LibmeshNodalShapeFunction.cpp
 * \author Stuart R. Slattery
 * \brief Nodal shape function implementation for STK mesh.
 */
//---------------------------------------------------------------------------//

#include "DTK_LibmeshNodalShapeFunction.hpp"
#include <DTK_DBC.hpp>

#include <libmesh/fe_compute_data.h>
#include <libmesh/fe_interface.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
LibmeshNodalShapeFunction::LibmeshNodalShapeFunction(
    const Teuchos::RCP<libMesh::MeshBase> &libmesh_mesh,
    const Teuchos::RCP<libMesh::System> &libmesh_system )
    : d_libmesh_mesh( libmesh_mesh )
    , d_libmesh_system( libmesh_system )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Given an entity, get the ids of its support locations
void LibmeshNodalShapeFunction::entitySupportIds(
    const Entity &entity, Teuchos::Array<SupportId> &support_ids ) const
{
    // Node case.
    if ( 0 == entity.topologicalDimension() )
    {
        DTK_CHECK( extractGeom<libMesh::Node>( entity )->valid_id() );
        support_ids.assign( 1, extractGeom<libMesh::Node>( entity )->id() );
    }

    // Element case.
    else
    {
        Teuchos::Ptr<libMesh::Elem> elem = extractGeom<libMesh::Elem>( entity );
        int num_nodes = elem->n_nodes();
        support_ids.resize( num_nodes );
        for ( int n = 0; n < num_nodes; ++n )
        {
            DTK_CHECK( elem->get_node( n )->valid_id() );
            support_ids[n] = elem->get_node( n )->id();
        }
    }
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the shape function of the
// entity at that point.
void LibmeshNodalShapeFunction::evaluateValue(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    Teuchos::Array<double> &values ) const
{
    int space_dim = entity.physicalDimension();
    libMesh::Point lm_reference_point;
    for ( int d = 0; d < space_dim; ++d )
    {
        lm_reference_point( d ) = reference_point[d];
    }

    libMesh::FEComputeData fe_compute_data(
        d_libmesh_system->get_equation_systems(), lm_reference_point );

    libMesh::FEInterface::compute_data(
        space_dim, d_libmesh_system->variable_type( 0 ),
        extractGeom<libMesh::Elem>( entity ).getRawPtr(), fe_compute_data );

    values = fe_compute_data.shape;
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the gradient of the shape
// function of the entity at that point.
void LibmeshNodalShapeFunction::evaluateGradient(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    Teuchos::Array<Teuchos::Array<double>> &gradients ) const
{
    return EntityShapeFunction::evaluateGradient( entity, reference_point,
                                                  gradients );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_LibmeshNodalShapeFunction.cpp
//---------------------------------------------------------------------------//
