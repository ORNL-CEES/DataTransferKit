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
 * \file DTK_BasicGeometryEntityImpl.hpp
 * \author Stuart R. Slattery
 * \brief BasicGeometryEntityImpl declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BASICGEOMETRYENTITYIMPL_HPP
#define DTK_BASICGEOMETRYENTITYIMPL_HPP

#include <iostream>

#include "DTK_EntityImpl.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class BasicGeometryEntityImpl
  \brief BasicGeometryEntityImpl interface.
  
  BasicGeometryEntityImpl gives an interface for simple geometries. These objects
  effectivelty define their own EntityImplLocalMap interface as these functions
  are typically statisfied with analytic expressions for basic geometric
  objects.
*/
//---------------------------------------------------------------------------//
class BasicGeometryEntityImpl : public EntityImpl
{
  public:

    // Default constructor.
    BasicGeometryEntityImpl();

    // Destructor.
    virtual ~BasicGeometryEntityImpl();

    //@{
    //! EntityImpl interface.
    /*!
     * \brief Get the entity type.
     * \return The entity type.
     */
    virtual EntityType entityType() const;

    /*!
     * \brief Get the unique global identifier for the entity.
     * \return A unique global identifier for the entity.
     */
    virtual EntityId id() const;
    
    /*!
     * \brief Get the parallel rank that owns the entity.
     * \return The parallel rank that owns the entity.
     */
    virtual int ownerRank() const;

    /*!
     * \brief Return the physical dimension of the entity.
     * \return The physical dimension of the entity. Any physical coordinates
     * describing the entity will be of this dimension.
     */
    virtual int physicalDimension() const;

    /*!
     * \brief Return the Cartesian bounding box around an entity.
     * \param bounds The bounds of the box
     * (x_min,y_min,z_min,x_max,y_max,z_max).
     */
    virtual void boundingBox( Teuchos::Tuple<double,6>& bounds ) const;

    /*!
     * \brief Determine if an entity is on the surface of the set.
     */
    virtual bool onSurface() const;

    /*!
     * \brief Determine if an entity is in the block with the given id.
     */
    virtual bool inBlock( const int block_id ) const;

    /*!
     * \brief Determine if an entity is on the boundary with the given id.
     */
    virtual bool onBoundary( const int boundary_id ) const;

    /*!
     * \brief Get the extra data on the entity.
     */
    virtual Teuchos::RCP<EntityExtraData> extraData() const;
    //@}

    //@{
    //! BasicGeometryEntityImpl interface.
    // Return the entity measure with respect to the parameteric
    virtual double measure() const;

    // Compute the centroid of the entity.
    virtual void centroid( const Teuchos::ArrayView<double>& centroid ) const;

    // (Safeguard the reverse map) Perform a safeguard check for mapping a
    // point to the reference space of an entity using the given tolerance.
    virtual bool isSafeToMapToReferenceFrame(
	const Teuchos::ArrayView<const double>& point ) const;

    // (Reverse Map) Map a point to the reference space of an entity. Return
    // the parameterized point.
    virtual bool mapToReferenceFrame( 
	const Teuchos::ArrayView<const double>& point,
	const Teuchos::ArrayView<double>& reference_point ) const;

    // Determine if a reference point is in the parameterized space of an
    // entity.
    virtual bool checkPointInclusion( 
	const double tolerance,
	const Teuchos::ArrayView<const double>& reference_point ) const;

    // (Forward Map) Map a reference point to the physical space of an entity.
    virtual void mapToPhysicalFrame( 
	const Teuchos::ArrayView<const double>& reference_point,
	const Teuchos::ArrayView<double>& point ) const;
    //@}
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_BASICGEOMETRYENTITYIMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_BasicGeometryEntityImpl.hpp
//---------------------------------------------------------------------------//

