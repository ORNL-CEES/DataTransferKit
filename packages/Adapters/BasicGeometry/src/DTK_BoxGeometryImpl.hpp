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
 * \file DTK_BoxGeometryImpl.hpp
 * \author Stuart R. Slattery
 * \brief BoxGeometry implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BOXGEOMETRYIMPL_HPP
#define DTK_BOXGEOMETRYIMPL_HPP

#include "DTK_BasicGeometryEntityImpl.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class BoxGeometryImpl
 * \brief Axis-aligned Cartesian box container implementation.
 *
 * All three dimensions are explictly represented in this bounding box. This
 * is different from a bounding box in that it must always be finite and of a
 * fixed 3 dimensions.
 */
//---------------------------------------------------------------------------//
class BoxGeometryImpl : public BasicGeometryEntityImpl
{
  public:
    // Default constructor.
    BoxGeometryImpl();

    // Constructor.
    BoxGeometryImpl( const EntityId global_id, const int owner_rank,
                     const int block_id, const double x_min, const double y_min,
                     const double z_min, const double x_max, const double y_max,
                     const double z_max );

    // Tuple constructor.
    BoxGeometryImpl( const EntityId global_id, const int owner_rank,
                     const int block_id,
                     const Teuchos::Tuple<double, 6> &bounds );

    // Get the unique global identifier for the entity.
    EntityId id() const override;

    // Get the parallel rank that owns the entity.
    int ownerRank() const override;

    // Return the topological dimension of the entity.
    int topologicalDimension() const override;

    // Return the physical dimension of the entity.
    int physicalDimension() const override;

    // Compute the bounding box around the entity.
    void boundingBox( Teuchos::Tuple<double, 6> &bounds ) const override;

    // Determine if an entity is in the block with the given id.
    bool inBlock( const int block_id ) const override;

    // Determine if an entity is on the boundary with the given id.
    bool onBoundary( const int boundary_id ) const override;

    // Provide a one line description of the object.
    std::string description() const override
    {
        return std::string( "Basic Geometry BoxGeometry" );
    }

    // Provide a verbose description of the object.
    void describe( Teuchos::FancyOStream &out,
                   const Teuchos::EVerbosityLevel verb_level ) const override;

    // Return the entity measure with respect to the parameteric
    double measure() const override;

    // Compute the centroid of the entity.
    void centroid( const Teuchos::ArrayView<double> &centroid ) const override;

    // (Reverse Map) Map a point to the reference space of an entity. Return
    // the parameterized point.
    bool mapToReferenceFrame(
        const Teuchos::ArrayView<const double> &point,
        const Teuchos::ArrayView<double> &reference_point ) const override;

    // Determine if a reference point is in the parameterized space of an
    // entity.
    bool checkPointInclusion( const double tolerance,
                              const Teuchos::ArrayView<const double>
                                  &reference_point ) const override;

    // (Forward Map) Map a reference point to the physical space of an entity.
    void mapToPhysicalFrame(
        const Teuchos::ArrayView<const double> &reference_point,
        const Teuchos::ArrayView<double> &point ) const override;

  private:
    // Global id.
    EntityId d_global_id;

    // Owning parallel rank.
    int d_owner_rank;

    // Block id.
    int d_block_id;

    // X min.
    double d_x_min;

    // Y min.
    double d_y_min;

    // Z min.
    double d_z_min;

    // X max.
    double d_x_max;

    // Y max.
    double d_y_max;

    // Z max.
    double d_z_max;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_BOXGEOMETRYIMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_BoxGeometry.hpp
//---------------------------------------------------------------------------//
