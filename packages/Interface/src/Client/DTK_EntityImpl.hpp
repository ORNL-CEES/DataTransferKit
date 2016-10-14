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
 * \brief DTK_EntityImpl.hpp
 * \author Stuart R. Slattery
 * \brief Geometric entity implementation.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYIMPL_HPP
#define DTK_ENTITYIMPL_HPP

#include "DTK_EntityExtraData.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class EntityImpl
  \brief Geometric entity implementation definition.
*/
//---------------------------------------------------------------------------//
class EntityImpl
{
  public:
    /*!
     * \brief Constructor.
     */
    EntityImpl() { /* ... */}

    /*!
     * \brief Destructor.
     */
    virtual ~EntityImpl() { /* ... */}

    /*!
     * \brief Get the unique global identifier for the entity.
     *
     * \return A unique global identifier for the entity.
     */
    virtual EntityId id() const = 0;

    /*!
     * \brief Get the parallel rank that owns the entity.
     *
     * \return The parallel rank that owns the entity.
     */
    virtual int ownerRank() const = 0;

    /*!
     * \brief Return the topological dimension of the entity.
     *
     * \return The topological dimension of the entity. Any parametric
     * coordinates describing the entity will be of this dimension.
     */
    virtual int topologicalDimension() const = 0;

    /*!
     * \brief Return the physical dimension of the entity.
     *
     * \return The physical dimension of the entity. Any physical coordinates
     * describing the entity will be of this dimension.
     */
    virtual int physicalDimension() const = 0;

    /*!
     * \brief Return the Cartesian bounding box around an entity.
     *
     * \param bounds The bounds of the box
     * (x_min,y_min,z_min,x_max,y_max,z_max).
     */
    virtual void boundingBox( Teuchos::Tuple<double, 6> &bounds ) const = 0;

    /*!
     * \brief Determine if an entity is in the block with the given id.
     */
    virtual bool inBlock( const int block_id ) const = 0;

    /*!
     * \brief Determine if an entity is on the boundary with the given id.
     */
    virtual bool onBoundary( const int boundary_id ) const = 0;

    /*!
     * \brief Get the extra data on the entity.
     */
    virtual Teuchos::RCP<EntityExtraData> extraData() const
    {
        return Teuchos::null;
    }

    /*!
     * \brief Provide a one line description of the object.
     */
    virtual std::string description() const
    {
        return std::string( "DataTransferKit::EntityImpl" );
    }

    /*!
     * \brief Provide a verbose description of the object.
     */
    virtual void describe( Teuchos::FancyOStream &out,
                           const Teuchos::EVerbosityLevel /*verb_level*/ ) const
    {
        out << "DataTransferKit::EntityImpl" << std::endl;
    }
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_ENTITYIMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_EntityImpl.hpp
//---------------------------------------------------------------------------//
