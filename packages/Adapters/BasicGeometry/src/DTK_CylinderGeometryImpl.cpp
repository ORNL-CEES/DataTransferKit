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
  * \file DTK_CylinderGeometryImpl.cpp
  * \author Stuart R. Slattery
  * \brief Bounding box definition.
  */
//---------------------------------------------------------------------------//

#include <cmath>

#include "DTK_CylinderGeometryImpl.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
CylinderGeometryImpl::CylinderGeometryImpl()
    : d_global_id( dtk_invalid_entity_id )
    , d_owner_rank( -1 )
    , d_block_id( 0 )
    , d_length( 0.0 )
    , d_radius( 0.0 )
    , d_centroid_x( 0.0 )
    , d_centroid_y( 0.0 )
    , d_centroid_z( 0.0 )
{ /* ... */
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param length Length of cylinder along Z-axis.
 *
 * \param radius Radius of cylinder.
 *
 * \param centroid_x Centroid X-coordinate.
 *
 * \param centroid_y Centroid Y-coordinate.
 *
 * \param centroid_z Centroid Z-coordinate.
 */
CylinderGeometryImpl::CylinderGeometryImpl(
    const EntityId global_id, const int owner_rank, const int block_id,
    const double length, const double radius, const double centroid_x,
    const double centroid_y, const double centroid_z )
    : d_global_id( global_id )
    , d_owner_rank( owner_rank )
    , d_block_id( block_id )
    , d_length( length )
    , d_radius( radius )
    , d_centroid_x( centroid_x )
    , d_centroid_y( centroid_y )
    , d_centroid_z( centroid_z )
{
    DTK_REQUIRE( 0.0 <= d_length );
    DTK_REQUIRE( 0.0 <= d_radius );
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
EntityId CylinderGeometryImpl::id() const { return d_global_id; }

//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
int CylinderGeometryImpl::ownerRank() const { return d_owner_rank; }

//---------------------------------------------------------------------------//
// Return the topological dimension of the entity.
int CylinderGeometryImpl::topologicalDimension() const { return 3; }

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
int CylinderGeometryImpl::physicalDimension() const { return 3; }

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the bounding box around the box.
 *
 * \return The bounding box bounds.
 */
void CylinderGeometryImpl::boundingBox(
    Teuchos::Tuple<double, 6> &bounds ) const
{
    bounds =
        Teuchos::tuple( d_centroid_x - d_radius, d_centroid_y - d_radius,
                        d_centroid_z - d_length / 2, d_centroid_x + d_radius,
                        d_centroid_y + d_radius, d_centroid_z + d_length / 2 );
}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
bool CylinderGeometryImpl::inBlock( const int block_id ) const
{
    return ( block_id == d_block_id );
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
bool CylinderGeometryImpl::onBoundary( const int boundary_id ) const
{
    return false;
}

//---------------------------------------------------------------------------//
void CylinderGeometryImpl::describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel /*verb_level*/ ) const
{
    out << "---" << std::endl;
    out << description() << std::endl;
    out << "Id: " << id() << std::endl;
    out << "Owner rank: " << ownerRank() << std::endl;
    out << "Block id: " << d_block_id << std::endl;
    out << "Length: " << d_length << std::endl;
    out << "Radius: " << d_radius << std::endl;
    out << "Centroid (x,y,z): " << d_centroid_x << " " << d_centroid_y << " "
        << d_centroid_z << std::endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the measure of the box.
 *
 * \return Return the measure of the box.
 */
double CylinderGeometryImpl::measure() const
{
    double zero = 0.0;
    double pi = 2.0 * std::acos( zero );
    return pi * d_radius * d_radius * d_length;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the centroid of the box.
 *
 * \return The centroid coordinates.
 */
void CylinderGeometryImpl::centroid(
    const Teuchos::ArrayView<double> &centroid ) const
{
    centroid[0] = d_centroid_x;
    centroid[1] = d_centroid_y;
    centroid[2] = d_centroid_z;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map a point to the reference space of an entity. Return the
 */
bool CylinderGeometryImpl::mapToReferenceFrame(
    const Teuchos::ArrayView<const double> &point,
    const Teuchos::ArrayView<double> &reference_point ) const
{
    reference_point.assign( point );
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine if a reference point is in the parameterized space of
 * an entity.
 */
bool CylinderGeometryImpl::checkPointInclusion(
    const double tolerance,
    const Teuchos::ArrayView<const double> &reference_point ) const
{
    DTK_REQUIRE( reference_point.size() == 3 );

    double x_dist = d_centroid_x - reference_point[0];
    double y_dist = d_centroid_y - reference_point[1];
    double r = std::sqrt( x_dist * x_dist + y_dist * y_dist );
    double rad_tol = d_radius * tolerance;
    double half_length_tol = d_length / 2 + d_length * tolerance;

    return ( ( r <= d_radius + rad_tol ) &&
             ( reference_point[2] >= d_centroid_z - half_length_tol ) &&
             ( reference_point[2] <= d_centroid_z + half_length_tol ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map a reference point to the physical space of an entity.
 */
void CylinderGeometryImpl::mapToPhysicalFrame(
    const Teuchos::ArrayView<const double> &reference_point,
    const Teuchos::ArrayView<double> &point ) const
{
    point.assign( reference_point );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_CylinderGeometryImpl.cpp
//---------------------------------------------------------------------------//
