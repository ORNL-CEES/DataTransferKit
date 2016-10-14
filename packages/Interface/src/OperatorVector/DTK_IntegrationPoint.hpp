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
 * \brief DTK_IntegrationPointSet.hpp
 * \author Stuart R. Slattery
 * \brief Integration point set.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTEGRATIONPOINT_HPP
#define DTK_INTEGRATIONPOINT_HPP

#include <limits>

#include "DTK_Entity.hpp"
#include "DTK_EntityImpl.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_Ptr.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// IntegrationPoint container.
//---------------------------------------------------------------------------//
class IntegrationPoint
{
  public:
    // Global id of the integration point.
    EntityId d_gid;

    // Measure of the entity that owns the point.
    double d_owner_measure;

    // Weight of the point in the integration rule.
    double d_integration_weight;

    // Physical coordinates of the point in the owning entity.
    Teuchos::Array<double> d_physical_coordinates;

    // Support ids of the owning entity.
    Teuchos::Array<SupportId> d_owner_support_ids;

    // Shape function evaluation of the point in the owning entity.
    Teuchos::Array<double> d_owner_shape_evals;
};

//---------------------------------------------------------------------------//
// IntegrationPoint EntityImpl subclass.
//---------------------------------------------------------------------------//
class IntegrationPointEntityImpl : public EntityImpl
{
  public:
    /*!
     * \brief Constructor
     */
    IntegrationPointEntityImpl( const Teuchos::Ptr<IntegrationPoint> &ip )
        : d_ip( ip )
    { /* ... */
    }

    /*!
     * \brief Get the unique global identifier for the entity.
     *
     * \return A unique global identifier for the entity.
     */
    EntityId id() const override { return d_ip->d_gid; }

    /*!
     * \brief Get the parallel rank that owns the entity.
     *
     * \return The parallel rank that owns the entity.
     */
    int ownerRank() const override { return -1; }

    /*!
     * \brief Return the topological dimension of the entity.
     *
     * \return The topological dimension of the entity. Any parametric
     * coordinates describing the entity will be of this dimension.
     */
    int topologicalDimension() const override { return 0; }

    /*!
     * \brief Return the physical dimension of the entity.
     *
     * \return The physical dimension of the entity. Any physical coordinates
     * describing the entity will be of this dimension.
     */
    int physicalDimension() const override
    {
        return d_ip->d_physical_coordinates.size();
    }

    /*!
     * \brief Return the Cartesian bounding box around an entity.
     *
     * \param bounds The bounds of the box
     * (x_min,y_min,z_min,x_max,y_max,z_max).
     */
    void boundingBox( Teuchos::Tuple<double, 6> &bounds ) const override
    {
        for ( int d = 0; d < physicalDimension(); ++d )
        {
            bounds[d] = d_ip->d_physical_coordinates[d];
            bounds[d + 3] = d_ip->d_physical_coordinates[d];
        }
        for ( int d = physicalDimension(); d < 3; ++d )
        {
            bounds[d] = -std::numeric_limits<double>::max();
            bounds[d + 3] = std::numeric_limits<double>::max();
        }
    }

    /*!
     * \brief Determine if an entity is in the block with the given id.
     */
    bool inBlock( const int block_id ) const override { return false; }

    /*!
     * \brief Determine if an entity is on the boundary with the given id.
     */
    bool onBoundary( const int boundary_id ) const override { return false; }

    /*!
     * \brief Provide a one line description of the object.
     */
    virtual std::string description() const override
    {
        return std::string( "DataTransferKit::IntegrationPointEntityImpl" );
    }

    /*!
     * \brief Provide a verbose description of the object.
     */
    virtual void
    describe( Teuchos::FancyOStream &out,
              const Teuchos::EVerbosityLevel /*verb_level*/ ) const override
    {
        out << this->description() << std::endl;
    }

  private:
    // Pointer to the underlying integration point.
    Teuchos::Ptr<IntegrationPoint> d_ip;
};

//---------------------------------------------------------------------------//
// IntegrationPoint Entity subclass.
//---------------------------------------------------------------------------//
class IntegrationPointEntity : public Entity
{
  public:
    IntegrationPointEntity( const Teuchos::Ptr<IntegrationPoint> &ip )
    {
        this->b_entity_impl =
            Teuchos::rcp( new IntegrationPointEntityImpl( ip ) );
    }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_INTEGRATIONPOINT_HPP

//---------------------------------------------------------------------------//
// end DTK_IntegrationPointSet.hpp
//---------------------------------------------------------------------------//
