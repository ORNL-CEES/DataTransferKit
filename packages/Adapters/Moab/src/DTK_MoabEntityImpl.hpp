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
 * \brief DTK_MoabEntityImpl.hpp
 * \author Stuart R. Slattery
 * \brief Moab entity implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOABENTITYIMPL_HPP
#define DTK_MOABENTITYIMPL_HPP

#include "DTK_Types.hpp"
#include "DTK_EntityImpl.hpp"
#include "DTK_MoabEntityExtraData.hpp"
#include "DTK_MoabMeshSetIndexer.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>

#include <MBParallelComm.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class MoabEntityImpl
  \brief Geometric entity implementation definition.
*/
//---------------------------------------------------------------------------//
class MoabEntityImpl : public EntityImpl
{
  public:

    /*!
     * \brief Constructor.
     */
    MoabEntityImpl( const moab::EntityHandle& moab_entity,
		    const Teuchos::Ptr<moab::ParallelComm>& moab_mesh,
		    const Teuchos::Ptr<MoabMeshSetIndexer>& set_indexer );

    /*!
     * \brief Destructor.
     */
    ~MoabEntityImpl();

    /*!
     * \brief Get the entity type.
     * \return The entity type.
     */
    EntityType entityType() const;

    /*!
     * \brief Get the unique global identifier for the entity.
     * \return A unique global identifier for the entity.
     */
    EntityId id() const;
    
    /*!
     * \brief Get the parallel rank that owns the entity.
     * \return The parallel rank that owns the entity.
     */
    int ownerRank() const;

    /*!
     * \brief Return the physical dimension of the entity.
     * \return The physical dimension of the entity. Any physical coordinates
     * describing the entity will be of this dimension.
     */
    int physicalDimension() const;

    /*!
     * \brief Return the Cartesian bounding box around an entity.
     * \param bounds The bounds of the box
     * (x_min,y_min,z_min,x_max,y_max,z_max).
     */
    void boundingBox( Teuchos::Tuple<double,6>& bounds ) const;

    /*!
     * \brief Determine if an entity is in the block with the given id.
     */
    bool inBlock( const int block_id ) const;

    /*!
     * \brief Determine if an entity is on the boundary with the given id.
     */
    bool onBoundary( const int boundary_id ) const;

    /*!
     * \brief Get the extra data on the entity.
     */
    Teuchos::RCP<EntityExtraData> extraData() const;

  private:

    // Moab entity extra data.
    Teuchos::RCP<MoabEntityExtraData> d_extra_data;

    // Moab parallel mesh.
    Teuchos::Ptr<moab::ParallelComm> d_moab_mesh;

    // Mesh set indexer.
    Teuchos::Ptr<MoabMeshSetIndexer> d_set_indexer;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MOABENTITYIMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_MoabEntityImpl.hpp
//---------------------------------------------------------------------------//
