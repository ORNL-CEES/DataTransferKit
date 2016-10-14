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
  * \file DTK_BoxGeometryImpl.cpp
  * \author Stuart R. Slattery
  * \brief BoxGeometry implementation.
  */
//---------------------------------------------------------------------------//

#include "DTK_BoxGeometryImpl.hpp"
#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Default constructor.
 */
BoxGeometryImpl::BoxGeometryImpl()
    : d_global_id( dtk_invalid_entity_id )
    , d_owner_rank( -1 )
    , d_block_id( 0 )
    , d_x_min( 0.0 )
    , d_y_min( 0.0 )
    , d_z_min( 0.0 )
    , d_x_max( 0.0 )
    , d_y_max( 0.0 )
    , d_z_max( 0.0 )
{ /* ... */
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param x_min Minimum x coordinate value in the box.
 *
 * \param y_min Minimum y coordinate value in the box.
 *
 * \param z_min Minimum z coordinate value in the box.
 *
 * \param x_max Maximum x coordinate value in the box.
 *
 * \param y_max Maximum y coordinate value in the box.
 *
 * \param z_max Maximum z coordinate value in the box.
 */
BoxGeometryImpl::BoxGeometryImpl( const EntityId global_id,
                                  const int owner_rank, const int block_id,
                                  const double x_min, const double y_min,
                                  const double z_min, const double x_max,
                                  const double y_max, const double z_max )
    : d_global_id( global_id )
    , d_owner_rank( owner_rank )
    , d_block_id( block_id )
    , d_x_min( x_min )
    , d_y_min( y_min )
    , d_z_min( z_min )
    , d_x_max( x_max )
    , d_y_max( y_max )
    , d_z_max( z_max )
{
    DTK_REQUIRE( d_x_min <= d_x_max );
    DTK_REQUIRE( d_y_min <= d_y_max );
    DTK_REQUIRE( d_z_min <= d_z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Tuple constructor.
 *
 * \param bounds Tuple containing {x_min, y_min, z_min, x_max, y_max, z_max}.
 */
BoxGeometryImpl::BoxGeometryImpl( const EntityId global_id,
                                  const int owner_rank, const int block_id,
                                  const Teuchos::Tuple<double, 6> &bounds )
    : d_global_id( global_id )
    , d_owner_rank( owner_rank )
    , d_block_id( block_id )
    , d_x_min( bounds[0] )
    , d_y_min( bounds[1] )
    , d_z_min( bounds[2] )
    , d_x_max( bounds[3] )
    , d_y_max( bounds[4] )
    , d_z_max( bounds[5] )
{
    DTK_REQUIRE( d_x_min <= d_x_max );
    DTK_REQUIRE( d_y_min <= d_y_max );
    DTK_REQUIRE( d_z_min <= d_z_max );
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
EntityId BoxGeometryImpl::id() const { return d_global_id; }

//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
int BoxGeometryImpl::ownerRank() const { return d_owner_rank; }

//---------------------------------------------------------------------------//
// Return the topological dimension of the entity.
int BoxGeometryImpl::topologicalDimension() const { return 3; }

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
int BoxGeometryImpl::physicalDimension() const { return 3; }

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the bounding box around the box.
 *
 * \return The bounding box bounds.
 */
void BoxGeometryImpl::boundingBox( Teuchos::Tuple<double, 6> &bounds ) const
{
    bounds =
        Teuchos::tuple( d_x_min, d_y_min, d_z_min, d_x_max, d_y_max, d_z_max );
}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
bool BoxGeometryImpl::inBlock( const int block_id ) const
{
    return ( block_id == d_block_id );
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
bool BoxGeometryImpl::onBoundary( const int boundary_id ) const
{
    return false;
}

//---------------------------------------------------------------------------//
void BoxGeometryImpl::describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel /*verb_level*/ ) const
{
    out << "---" << std::endl;
    out << description() << std::endl;
    out << "Id: " << id() << std::endl;
    out << "Owner rank: " << ownerRank() << std::endl;
    out << "Block id: " << d_block_id << std::endl;
    out << "(xmin,ymin,zmin,xmax,ymax,zmax): " << d_x_min << " " << d_y_min
        << " " << d_z_min << " " << d_x_max << " " << d_y_max << " " << d_z_max
        << std::endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the measure of the box.
 *
 * \return Return the measure of the box.
 */
double BoxGeometryImpl::measure() const
{
    return ( d_x_max - d_x_min ) * ( d_y_max - d_y_min ) *
           ( d_z_max - d_z_min );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the centroid of the box.
 *
 * \return The centroid coordinates.
 */
void BoxGeometryImpl::centroid(
    const Teuchos::ArrayView<double> &centroid ) const
{
    centroid[0] = ( d_x_max + d_x_min ) / 2.0;
    centroid[1] = ( d_y_max + d_y_min ) / 2.0;
    centroid[2] = ( d_z_max + d_z_min ) / 2.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map a point to the reference space of an entity. Return the
 */
bool BoxGeometryImpl::mapToReferenceFrame(
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
bool BoxGeometryImpl::checkPointInclusion(
    const double tolerance,
    const Teuchos::ArrayView<const double> &reference_point ) const
{
    DTK_REQUIRE( 3 == reference_point.size() );

    double x_tol = ( d_x_max - d_x_min ) * tolerance;
    double y_tol = ( d_y_max - d_y_min ) * tolerance;
    double z_tol = ( d_z_max - d_z_min ) * tolerance;

    if ( reference_point[0] >= d_x_min - x_tol &&
         reference_point[1] >= d_y_min - y_tol &&
         reference_point[2] >= d_z_min - z_tol &&
         reference_point[0] <= d_x_max + x_tol &&
         reference_point[1] <= d_y_max + y_tol &&
         reference_point[2] <= d_z_max + z_tol )
    {
        return true;
    }

    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map a reference point to the physical space of an entity.
 */
void BoxGeometryImpl::mapToPhysicalFrame(
    const Teuchos::ArrayView<const double> &reference_point,
    const Teuchos::ArrayView<double> &point ) const
{
    point.assign( reference_point );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_BoxGeometryImpl.cpp
//---------------------------------------------------------------------------//
