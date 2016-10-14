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
 * \file DTK_CoarseGlobalSearch.hpp
 * \author Stuart R. Slattery
 * \brief CoarseGlobalSearch declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_COARSEGLOBALSEARCH_HPP
#define DTK_COARSEGLOBALSEARCH_HPP

#include "DTK_DBC.hpp"
#include "DTK_EntityIterator.hpp"
#include "DTK_EntityLocalMap.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class CoarseGlobalSearch
 * \brief A CoarseGlobalSearch data structure for global entity coarse search.
 */
//---------------------------------------------------------------------------//
class CoarseGlobalSearch
{
  public:
    // Constructor.
    CoarseGlobalSearch( const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                        const int physical_dimension,
                        const EntityIterator &domain_iterator,
                        const Teuchos::ParameterList &parameters );

    // Redistribute a set of range entity centroid coordinates with their
    // owner ranks to the owning domain process.
    void search( const EntityIterator &range_iterator,
                 const Teuchos::RCP<EntityLocalMap> &range_local_map,
                 const Teuchos::ParameterList &parameters,
                 Teuchos::Array<EntityId> &range_entity_ids,
                 Teuchos::Array<int> &range_owner_ranks,
                 Teuchos::Array<double> &range_centroids ) const;

    /*!
     * \brief Return the ids of the range entities that were not during the
     * last search (i.e. those that are guaranteed to not receive data from
     * the transfer).
     *
     * \return A view of the ids.
     */
    Teuchos::ArrayView<const EntityId> getMissedRangeEntityIds() const;

  private:
    // Assemble the local bounding box around an iterator.
    void assembleBoundingBox( const EntityIterator &entity_iterator,
                              Teuchos::Tuple<double, 6> &bounding_box ) const;

    // Check if two bounding boxes have an intersection.
    inline bool boxesIntersect( const Teuchos::Tuple<double, 6> &box_A,
                                const Teuchos::Tuple<double, 6> &box_B,
                                const double tolerance ) const;

    // Determine if a point is in a bounding box.
    inline bool pointInBox( const Teuchos::ArrayView<const double> &point,
                            const Teuchos::Tuple<double, 6> &box,
                            const double tolerance ) const;

  private:
    // Communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> d_comm;

    // Spatial dimension.
    int d_space_dim;

    // Domain bounding boxes.
    Teuchos::Array<Teuchos::Tuple<double, 6>> d_domain_boxes;

    // Boolean for tracking missed range entities.
    bool d_track_missed_range_entities;

    // An array of range entity ids that were not mapped during the last call
    // to setup.
    mutable Teuchos::Array<EntityId> d_missed_range_entity_ids;

    // Point inclusion tolerance.
    double d_inclusion_tol;
};

//---------------------------------------------------------------------------//
// Inline functions.
//---------------------------------------------------------------------------//
// Check if two bounding boxes have an intersection.
bool CoarseGlobalSearch::boxesIntersect( const Teuchos::Tuple<double, 6> &box_A,
                                         const Teuchos::Tuple<double, 6> &box_B,
                                         const double tolerance ) const
{
    double x_tol_A = ( box_A[3] - box_A[0] ) * tolerance;
    double y_tol_A = ( box_A[4] - box_A[1] ) * tolerance;
    double z_tol_A = ( box_A[5] - box_A[2] ) * tolerance;

    double x_tol_B = ( box_B[3] - box_B[0] ) * tolerance;
    double y_tol_B = ( box_B[4] - box_B[1] ) * tolerance;
    double z_tol_B = ( box_B[5] - box_B[2] ) * tolerance;

    return !( ( ( box_A[0] - x_tol_A ) > ( box_B[3] + x_tol_B ) ||
                ( box_A[3] + x_tol_A ) < ( box_B[0] - x_tol_B ) ) ||
              ( ( box_A[1] - y_tol_A ) > ( box_B[4] + y_tol_B ) ||
                ( box_A[4] + y_tol_A ) < ( box_B[1] - y_tol_B ) ) ||
              ( ( box_A[2] - z_tol_A ) > ( box_B[5] + z_tol_B ) ||
                ( box_A[5] + z_tol_A ) < ( box_B[2] - z_tol_B ) ) );
}

//---------------------------------------------------------------------------//
// Determine if a point is in a bounding box.
bool CoarseGlobalSearch::pointInBox(
    const Teuchos::ArrayView<const double> &point,
    const Teuchos::Tuple<double, 6> &box, const double tolerance ) const
{
    double x_tol = ( box[3] - box[0] ) * tolerance;
    double y_tol = ( box[4] - box[1] ) * tolerance;
    double z_tol = ( box[5] - box[2] ) * tolerance;

    if ( 3 == point.size() )
    {
        if ( point[0] >= ( box[0] - x_tol ) && point[1] >= ( box[1] - y_tol ) &&
             point[2] >= ( box[2] - z_tol ) && point[0] <= ( box[3] + x_tol ) &&
             point[1] <= ( box[4] + y_tol ) && point[2] <= ( box[5] + z_tol ) )
        {
            return true;
        }
    }
    else if ( 2 == point.size() )
    {
        if ( point[0] >= ( box[0] - x_tol ) && point[1] >= ( box[1] - y_tol ) &&
             point[0] <= ( box[3] + x_tol ) && point[1] <= ( box[4] + y_tol ) )
        {
            return true;
        }
    }
    else if ( 1 == point.size() )
    {
        if ( point[0] >= ( box[0] - x_tol ) && point[0] <= ( box[3] + x_tol ) )
        {
            return true;
        }
    }

    return false;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // DTK_COARSEGLOBALSEARCH_HPP

//---------------------------------------------------------------------------//
// end CoarseGlobalSearch.hpp
//---------------------------------------------------------------------------//
