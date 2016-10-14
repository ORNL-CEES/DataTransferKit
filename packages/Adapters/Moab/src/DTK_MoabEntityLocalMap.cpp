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
 * \brief DTK_MoabEntityLocalMap.cpp
 * \author Stuart R. Slattery
 * \brief Forward and reverse local mappings for entities.
 */
//---------------------------------------------------------------------------//

#include "DTK_MoabEntityLocalMap.hpp"
#include "DTK_DBC.hpp"
#include "DTK_MoabHelpers.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
MoabEntityLocalMap::MoabEntityLocalMap(
    const Teuchos::RCP<moab::ParallelComm> &moab_mesh )
    : d_moab_mesh( moab_mesh )
    , d_inclusion_tol( 1.0e-6 )
    , d_newton_tol( 1.0e-9 )
{
    d_moab_evaluator =
        Teuchos::rcp( new moab::ElemEvaluator( d_moab_mesh->get_moab() ) );
}

//---------------------------------------------------------------------------//
// Set parameters for mapping.
void MoabEntityLocalMap::setParameters(
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
double MoabEntityLocalMap::measure( const Entity &entity ) const
{
    cacheEntity( entity );

    Teuchos::Array<double> measure( 3, 0.0 );
    DTK_CHECK_ERROR_CODE( d_moab_evaluator->integrate( measure.getRawPtr() ) );

    return measure[0];
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
void MoabEntityLocalMap::centroid(
    const Entity &entity, const Teuchos::ArrayView<double> &centroid ) const
{
    // Node case.
    if ( 0 == entity.topologicalDimension() )
    {
        moab::EntityHandle handle = MoabHelpers::extractEntity( entity );
        d_moab_mesh->get_moab()->get_coords( &handle, 1, centroid.getRawPtr() );
    }
    // Element case.
    else
    {
        cacheEntity( entity );
        Teuchos::Array<double> param_center;
        parametricCenter( entity, param_center );

        DTK_CHECK_ERROR_CODE( d_moab_evaluator->eval( param_center.getRawPtr(),
                                                      centroid.getRawPtr() ) );
    }
}

//---------------------------------------------------------------------------//
// Perform a safeguard check for mapping a point to the reference space
// of an entity using the given tolerance.
bool MoabEntityLocalMap::isSafeToMapToReferenceFrame(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &physical_point ) const
{
    int space_dim = entity.physicalDimension();
    int param_dim = d_moab_mesh->get_moab()->dimension_from_handle(
        MoabHelpers::extractEntity( entity ) );
    if ( space_dim == param_dim )
    {
        return EntityLocalMap::isSafeToMapToReferenceFrame( entity,
                                                            physical_point );
    }
    else
    {
        bool not_implemented = true;
        DTK_INSIST( !not_implemented );
    }
    return false;
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the parameterized
// point.
bool MoabEntityLocalMap::mapToReferenceFrame(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &physical_point,
    const Teuchos::ArrayView<double> &reference_point ) const
{
    cacheEntity( entity );

    int is_inside = -1;

    // Ignore the error code on this one because of the ridiculous
    // tolerancing/convergence scheme used in the implementation.
    d_moab_evaluator->reverse_eval( physical_point.getRawPtr(), d_newton_tol,
                                    d_inclusion_tol,
                                    reference_point.getRawPtr(), &is_inside );
    return ( is_inside > 0 );
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of an entity.
bool MoabEntityLocalMap::checkPointInclusion(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point ) const
{
    cacheEntity( entity );

    int is_inside = d_moab_evaluator->inside( reference_point.getRawPtr(),
                                              d_inclusion_tol );
    return ( is_inside > 0 );
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
void MoabEntityLocalMap::mapToPhysicalFrame(
    const Entity &entity,
    const Teuchos::ArrayView<const double> &reference_point,
    const Teuchos::ArrayView<double> &physical_point ) const
{
    cacheEntity( entity );

    DTK_CHECK_ERROR_CODE( d_moab_evaluator->eval(
        reference_point.getRawPtr(), physical_point.getRawPtr() ) );
}

//---------------------------------------------------------------------------//
// Compute the normal on a face (3D) or edge (2D) at a given reference point.
void MoabEntityLocalMap::normalAtReferencePoint(
    const Entity &entity, const Entity &parent_entity,
    const Teuchos::ArrayView<const double> &reference_point,
    const Teuchos::ArrayView<double> &normal ) const
{
    bool not_implemented = true;
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Cache an entity in the evaluator.
void MoabEntityLocalMap::cacheEntity( const Entity &entity ) const
{
    DTK_CHECK_ERROR_CODE( d_moab_evaluator->set_eval_set(
        MoabHelpers::extractEntity( entity ) ) );
    DTK_CHECK_ERROR_CODE( d_moab_evaluator->set_ent_handle(
        MoabHelpers::extractEntity( entity ) ) );
    DTK_CHECK_ERROR_CODE( d_moab_evaluator->set_tag( "COORDS", 0 ) );
}

//---------------------------------------------------------------------------//
// Get the parameteric center of an entity.
void MoabEntityLocalMap::parametricCenter(
    const Entity &entity, Teuchos::Array<double> &center ) const
{
    moab::EntityType moab_type = d_moab_mesh->get_moab()->type_from_handle(
        MoabHelpers::extractEntity( entity ) );

    switch ( moab_type )
    {
    case moab::MBTRI:
        center.assign( 2, 1.0 / 3.0 );
        break;
    case moab::MBQUAD:
        center.assign( 2, 0.0 );
        break;
    case moab::MBTET:
        center.assign( 3, 1.0 / 6.0 );
        break;
    case moab::MBHEX:
        center.assign( 3, 0.0 );
        break;
    default:
        center.resize( 0 );
        break;
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_MoabEntityLocalMap.cpp
//---------------------------------------------------------------------------//
