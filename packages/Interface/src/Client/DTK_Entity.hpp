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
 * \brief DTK_Entity.hpp
 * \author Stuart R. Slattery
 * \brief Geometric entity interface.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITY_HPP
#define DTK_ENTITY_HPP

#include "DTK_EntityExtraData.hpp"
#include "DTK_EntityImpl.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class Entity
  \brief Geometric entity interface definition.
*/
//---------------------------------------------------------------------------//
class Entity : public Teuchos::Describable
{
  public:
    /*!
     * \brief Constructor.
     */
    Entity();

    /*!
     * \brief Copy constructor.
     */
    Entity( const Entity &rhs );

    /*!
     *\brief Copy assignment operator.
     */
    Entity &operator=( const Entity &rhs );

    /*!
     * \brief Move constructor.
     */
    Entity( Entity &&rhs );

    /*!
     *\brief Move assignment operator.
     */
    Entity &operator=( Entity &&rhs );

    /*!
     * \brief Destructor.
     */
    virtual ~Entity();

    //@{
    //! Client interface.
    /*!
     * \brief Get the unique global identifier for the entity. Entities of
     * different topological dimensions may have the same id.
     *
     * \return A unique global identifier for the entity.
     */
    EntityId id() const;

    /*!
     * \brief Get the parallel rank that owns the entity.
     *
     * \return The parallel rank that owns the entity.
     */
    int ownerRank() const;

    /*!
     * \brief Return the topological dimension of the entity.
     *
     * \return The topological dimension of the entity. This is the dimension
     * of the entity reference frame. Any parametric coordinates describing
     * the entity will be of this dimension.
     */
    int topologicalDimension() const;

    /*!
     * \brief Return the physical dimension of the entity.
     *
     * \return The physical dimension of the entity. This is the dimension of
     * the entity physical frame. Any physical coordinates describing the
     * entity will be of this dimension.
     */
    int physicalDimension() const;

    /*!
     * \brief Return the Cartesian bounding box around an entity.
     *
     * \param bounds The bounds of the box
     * (x_min,y_min,z_min,x_max,y_max,z_max).
     */
    void boundingBox( Teuchos::Tuple<double, 6> &bounds ) const;

    /*!
     * \brief Determine if an entity is in the block with the given id.
     */
    bool inBlock( const int block_id ) const;

    /*!
     * \brief Determine if an entity is on the boundary with the given id.
     */
    bool onBoundary( const int boundary_id ) const;

    /*!
     * \brief Get the extra data on the entity. This is a convenient helper
     * for implementing the other interfaces.
     */
    Teuchos::RCP<EntityExtraData> extraData() const;
    //@}

    //@{
    //! Teuchos::Describable interface.
    /*!
     * \brief Provide a one line description of the object.
     */
    std::string description() const override;

    /*!
     * \brief Provide a verbose description of the object.
     */
    void describe( Teuchos::FancyOStream &out,
                   const Teuchos::EVerbosityLevel verb_level =
                       Teuchos::Describable::verbLevel_default ) const override;
    //@}

  protected:
    // Entity implementation.
    Teuchos::RCP<EntityImpl> b_entity_impl;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_ENTITY_HPP

//---------------------------------------------------------------------------//
// end DTK_Entity.hpp
//---------------------------------------------------------------------------//
