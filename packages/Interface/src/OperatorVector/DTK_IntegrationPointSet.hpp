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

#ifndef DTK_INTEGRATIONPOINTSET_HPP
#define DTK_INTEGRATIONPOINTSET_HPP

#include "DTK_EntityIterator.hpp"
#include "DTK_EntityLocalMap.hpp"
#include "DTK_IntegrationPoint.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class IntegrationPointSetIterator
  \brief ementation of iterator over entities in a basic set.
*/
class IntegrationPointSetIterator : public EntityIterator
{
  public:
    // Default constructor.
    IntegrationPointSetIterator();

    // Constructor.
    IntegrationPointSetIterator(
        Teuchos::RCP<Teuchos::Array<IntegrationPoint>> points );

    // Copy constructor.
    IntegrationPointSetIterator( const IntegrationPointSetIterator &rhs );

    /*!
     * \brief Assignment operator.
     */
    IntegrationPointSetIterator &
    operator=( const IntegrationPointSetIterator &rhs );

    // Pre-increment operator.
    EntityIterator &operator++() override;

    // Dereference operator.
    Entity &operator*(void)override;

    // Dereference operator.
    Entity *operator->(void)override;

    // Equal comparison operator.
    bool operator==( const EntityIterator &rhs ) const override;

    // Not equal comparison operator.
    bool operator!=( const EntityIterator &rhs ) const override;

    // An iterator assigned to the beginning.
    EntityIterator begin() const override;

    // An iterator assigned to the end.
    EntityIterator end() const override;

    // Create a clone of the iterator. We need this for the copy constructor
    // and assignment operator to pass along the underlying implementation.
    std::unique_ptr<EntityIterator> clone() const override;

  private:
    // Map to iterate over.
    Teuchos::RCP<Teuchos::Array<IntegrationPoint>> d_points;

    // Iterator over the entity map.
    Teuchos::Array<IntegrationPoint>::iterator d_points_it;

    // The current entity.
    Entity d_current_entity;
};

//---------------------------------------------------------------------------//
/*!
  \class IntegrationPointSet
  \brief EntitySet of integration points.
*/
//---------------------------------------------------------------------------//
class IntegrationPointSet : public EntityLocalMap
{
  public:
    /*!
     * \brief Constructor.
     */
    IntegrationPointSet( const Teuchos::RCP<const Teuchos::Comm<int>> &comm );

    // Add an integration point to the set.
    void addPoint( const IntegrationPoint &ip );

    // Finalize the point set to construct global ids.
    void finalize();

    // Get an integration point with the given global id.
    const IntegrationPoint &getPoint( const EntityId ip_id ) const;

    // Get an entity iterator over the integration points.
    EntityIterator entityIterator() const;

    // Get the number of points.
    int numPoints() const { return d_points.size(); }

    // Get the global maximum support size for all integration points.
    int globalMaxSupportSize() const;

    //@{
    //! EntityLocalMap interface. Only implement the centroid function for the
    // search.
    void centroid( const Entity &entity,
                   const Teuchos::ArrayView<double> &centroid ) const override;

    //! Not implemented for IntegrationPoints
    void setParameters( const Teuchos::ParameterList &parameters ) override
    { /* ...*/}
    double measure( const Entity &entity ) const override { return 0.0; }
    bool mapToReferenceFrame(
        const Entity &entity,
        const Teuchos::ArrayView<const double> &physical_point,
        const Teuchos::ArrayView<double> &reference_point ) const override
    {
        return false;
    }
    bool checkPointInclusion(
        const Entity &entity,
        const Teuchos::ArrayView<const double> &reference_point ) const override
    {
        return false;
    }
    void mapToPhysicalFrame(
        const Entity &entity,
        const Teuchos::ArrayView<const double> &reference_point,
        const Teuchos::ArrayView<double> &physical_point ) const override
    { /* ...*/
    }
    //@}

  private:
    // Communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> d_comm;

    // Local integration points.
    mutable Teuchos::Array<IntegrationPoint> d_points;

    // Starting global id for this proc.
    EntityId d_start_gid;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_INTEGRATIONPOINTSET_HPP

//---------------------------------------------------------------------------//
// end DTK_IntegrationPointSet.hpp
//---------------------------------------------------------------------------//
