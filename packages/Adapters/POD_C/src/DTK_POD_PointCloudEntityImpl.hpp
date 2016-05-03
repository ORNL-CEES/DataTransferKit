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
 * \brief DTK_POD_PointCloudEntityImpl.hpp
 * \author Stuart R. Slattery
 * \brief Point cloud entity implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_POD_POINTCLOUDENTITYIMPL_HPP
#define DTK_POD_POINTCLOUDENTITYIMPL_HPP

#include "DTK_Types.hpp"
#include "DTK_EntityImpl.hpp"
#include "DTK_C_API.h"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class POD_PointCloudEntityImpl
  \brief Geometric entity implementation definition.
*/
//---------------------------------------------------------------------------//
class POD_PointCloudEntityImpl : public EntityImpl
{
  public:

    /*!
     * \brief Constructor.
     */
    POD_PointCloudEntityImpl( const double* cloud_coords,
                              const unsigned num_points,
                              const int space_dim,
                              const DTK_Data_layout layout,
                              const EntityId global_id,                          
                              const int local_id,
                              const int owner_rank );

    /*
     * \brief Get the coordinates of the point in a given dimension.
     */
    double coords( const int dim ) const;

    //@{
    //! EntityImpl Interface.
    /*!
     * \brief Get the unique global identifier for the entity.
     * \return A unique global identifier for the entity.
     */
    EntityId id() const override;
    
    /*!
     * \brief Get the parallel rank that owns the entity.
     * \return The parallel rank that owns the entity.
     */
    int ownerRank() const override;

    /*!
     * \brief Return the topological dimension of the entity.  
     *
     * \return The topological dimension of the entity. Any parametric
     * coordinates describing the entity will be of this dimension.
     */
    int topologicalDimension() const override;

        /*!
     * \brief Return the physical dimension of the entity.
     * \return The physical dimension of the entity. Any physical coordinates
     * describing the entity will be of this dimension.
     */
    int physicalDimension() const override;

    /*!
     * \brief Return the Cartesian bounding box around an entity.
     * \param bounds The bounds of the box
     * (x_min,y_min,z_min,x_max,y_max,z_max).
     */
    void boundingBox( Teuchos::Tuple<double,6>& bounds ) const override;

    /*!
     * \brief Determine if an entity is in the block with the given id.
     */
    bool inBlock( const int block_id ) const override;

    /*!
     * \brief Determine if an entity is on the boundary with the given id.
     */
    bool onBoundary( const int boundary_id ) const override;

    /*!
     * \brief Get the extra data on the entity.
     */
    Teuchos::RCP<EntityExtraData> extraData() const override;

    /*!
     * \brief Provide a one line description of the object.
     */
    std::string description() const override
    { return std::string("POD Point Cloud Entity"); }

    /*!
     * \brief Provide a verbose description of the object.
     */
    void describe(
	Teuchos::FancyOStream& out,
	const Teuchos::EVerbosityLevel verb_level ) const override;
    //@}

  private:

    // Point cloud coordinates
    const double* d_cloud_coords;

    // Coordinate offsets.
    Teuchos::Array<int> d_offsets;

    // Global id.
    EntityId d_global_id;

    // Owner rank.
    int d_owner_rank;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_POD_POINTCLOUDENTITYIMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_POD_PointCloudEntityImpl.hpp
//---------------------------------------------------------------------------//
