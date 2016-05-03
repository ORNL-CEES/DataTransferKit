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
 * \brief DTK_POD_PointCloudEntitySet.hpp
 * \author Stuart R. Slattery
 * \brief POD point cloud entity set.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_POD_POINTCLOUDENTITYSET_HPP
#define DTK_POD_POINTCLOUDENTITYSET_HPP

#include <functional>

#include "DTK_EntitySet.hpp"
#include "DTK_Types.hpp"
#include "DTK_Entity.hpp"
#include "DTK_EntityIterator.hpp"

#include "DTK_POD_Types.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class POD_PointCloudEntitySet
  \brief Entity set implementation for point cloud.
*/
//---------------------------------------------------------------------------//
class POD_PointCloudEntitySet : public EntitySet
{
  public:

    /*!
     * \brief Constructor.
     */
    POD_PointCloudEntitySet(
        const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
        const double* cloud_coords,
        const EntityId* global_ids,
        const unsigned num_points,
        const int space_dim,
        const DataLayout layout );

    /*!
     * \brief Get the parallel communicator for the entity set.
     * \return A reference-counted pointer to the parallel communicator.
     */
    Teuchos::RCP<const Teuchos::Comm<int> > communicator() const override;
    //@}

    //@{
    //! Geometric data functions.
    /*!
     * \brief Return the largest physical dimension of the entities in the
     * set.  
     * \return The physical dimension of the set.
     */
    int physicalDimension() const override;
    //@}

    //@{
    //! Entity access functions.
    /*!
     * \brief Given an EntityId, get the entity.
     * \param entity_id Get the entity with this id.
     * \param topological_dimension Get the entity with this topological
     * dimension.
     * \param entity The entity with the given id.
     */
    void getEntity( const EntityId entity_id,
		    const int topological_dimension,
		    Entity& entity ) const override;

    /*!
     * \brief Get a iterator of the given entity type that satisfy the given
     * predicate.
     * \param topological_dimension The topological dimension of entity to get
     * an iterator for.
     * \param predicate The selection predicate.
     * \return A iterator of entities of the given type.
     */
    EntityIterator entityIterator( 
	const int topological_dimension,
	const PredicateFunction& predicate ) const override;

    /*!
     * \brief Given an entity, get the entities of the given topological
     * dimension that are adjacent to it.
     */
    void getAdjacentEntities(
	const Entity& entity,
	const int adjacent_dimension,
	Teuchos::Array<Entity>& adjacent_entities ) const override;

  private:

    // Communicator.
    Teuchos::RCP<const Teuchos::Comm<int> > d_comm;

    // Point cloud coordinates.
    const double* d_cloud_coords;

    // Point global ids.
    const EntityId* d_global_ids;
    
    // Number of points in the point cloud.
    unsigned d_num_points;

    // Spatial dimension of point cloud.
    int d_space_dim;

    // Layout of the point cloud.
    DataLayout d_layout;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_POD_POINTCLOUDENTITYSET_HPP

//---------------------------------------------------------------------------//
// end DTK_POD_PointCloudEntitySet.hpp
//---------------------------------------------------------------------------//
