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
 * \brief DTK_EntityLocalMap.hpp
 * \author Stuart R. Slattery
 * \brief Forward and reverse local mappings for entities.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_ENTITYLOCALMAP_HPP
#define DTK_ENTITYLOCALMAP_HPP

#include "DTK_Entity.hpp"
#include "DTK_Types.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class EntityLocalMap
  \brief Entity forward and reverse local map interface definition.

  An EntityLocalMap provides an interface for accessing forward and reverse
  maps for an entity's local coordinates as well as related convenience
  functions.
*/
//---------------------------------------------------------------------------//
class EntityLocalMap
{
  public:
    /*!
     * \brief Constructor.
     */
    EntityLocalMap();

    /*!
     * \brief Destructor.
     */
    virtual ~EntityLocalMap();

    /*
     * \brief Set parameters for mapping.
     *
     * \param parameters Parameters for mapping.
     */
    virtual void setParameters( const Teuchos::ParameterList &parameters ) = 0;

    /*!
     * \brief Return the entity measure in the physical frame with respect to
     * the parameteric dimension (volume for a 3D entity, area for 2D, length
     * for 1D, 0 for 0D).
     *
     * \param entity Compute the measure for this entity.
     *
     * \return The measure of the entity.
     */
    virtual double measure( const Entity &entity ) const = 0;

    /*!
     * \brief Return the centroid of the entity in the physical frame.
     *
     * \param centroid A view of the centroid coordinates. This view will
     * be allocated. Assign a view of your centroid to this view.
     */
    virtual void
    centroid( const Entity &entity,
              const Teuchos::ArrayView<double> &centroid ) const = 0;

    /*!
     * \brief (Safeguard the reverse map) Perform a safeguard check for
     * mapping a point to the reference space of an entity using the given
     * tolerance.
     *
     * Mapping from physical domain to a reference one requires solving a
     * nonlinear problem. It has been observed (see [1], sections VI.2.1-VI.2.3)
     * that this nonlinear method may not converge for some points that should
     * not have been considered in the first place as they are too far. The
     * method tries to address the issue by providing a way to ignore those
     * points early in the process. The default implementation checks if the
     * physical point is in the bounding box of the entity.
     *
     * It is recommended that a user calls isSafeToMapToReferenceFrame() before
     * calling mapToReferenceFrame(). The tolerance is provided by a specific
     * implementation of the interface.
     *
     * [1] Lebrun-Grandie, Damien Thomas (2014). Contact Detection and
     * Constraints Enforcement for the Simulation of Pellet/Clad
     * Thermo-Mechanical Contact in Nuclear Fuel Rods. Doctoral dissertation,
     * Texas A & M University
     *
     * \param entity Perfrom the mapping for this entity.
     *
     * \param parameters Parameters to be used for the safeguard check.
     *
     * \param physical_point A view into an array of size physicalDimension()
     * containing the coordinates of the point to map.
     *
     * \return Return true if it is safe to map to the reference frame.
     */
    virtual bool isSafeToMapToReferenceFrame(
        const Entity &entity,
        const Teuchos::ArrayView<const double> &physical_point ) const;

    /*!
     * \brief (Reverse Map) Map a point to the reference space of an
     * entity. Return the parameterized point.
     *
     * \param entity Perfrom the mapping for this entity.
     *
     * \param parameters Parameters to be used for the mapping procedure.
     *
     * \param physical_point A view into an array of size physicalDimension()
     * containing the coordinates of the point to map.
     *
     * \param reference_point A view into an array of size physicalDimension()
     * to write the reference coordinates of the mapped point.
     *
     * \return Return true if the map to reference frame succeeded.
     */
    virtual bool mapToReferenceFrame(
        const Entity &entity,
        const Teuchos::ArrayView<const double> &physical_point,
        const Teuchos::ArrayView<double> &reference_point ) const = 0;

    /*!
     * \brief Determine if a reference point is in the parameterized space of
     * an entity.
     *
     * \param entity Perfrom the mapping for this entity.
     *
     * \param parameters Parameters to be used for the point inclusion check.
     *
     * \param reference_point A view into an array of size physicalDimension()
     * containing the reference coordinates of the mapped point.
     *
     * \return True if the point is in the reference space, false if not.
     */
    virtual bool checkPointInclusion(
        const Entity &entity,
        const Teuchos::ArrayView<const double> &reference_point ) const = 0;

    /*!
     * \brief (Forward Map) Map a reference point to the physical space of an
     * entity.
     *
     * \param entity Perfrom the mapping for this entity.
     *
     * \param reference_point A view into an array of size physicalDimension()
     * containing the reference coordinates of the mapped point.
     *
     * \param physical_point A view into an array of size physicalDimension()
     * to write the coordinates of physical point.
     */
    virtual void mapToPhysicalFrame(
        const Entity &entity,
        const Teuchos::ArrayView<const double> &reference_point,
        const Teuchos::ArrayView<double> &physical_point ) const = 0;

    /*!
     * \brief Compute the normal on a face (3D) or edge (2D) at a given
     * reference point. A default implementation is provided using a finite
     * difference scheme.
     *
     * \param entity Compute the normal for this entity.
     *
     * \param parent_entity The adjacent parent entity used to determine which
     * direction is outward. The parent entity should be of a higher
     * topological dimension than the entity and be adjacent to the entity.
     *
     * \param reference_point Compute the normal at this reference point.
     *
     * \param normal A view into an array of size physicalDimension() to write
     * the normal.
     *
     * \throw DataTransferKitException The function throws if a normal cannot be
     * calculated
     */
    virtual void normalAtReferencePoint(
        const Entity &entity, const Entity &parent_entity,
        const Teuchos::ArrayView<const double> &reference_point,
        const Teuchos::ArrayView<double> &normal ) const;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_ENTITYLOCALMAP_HPP

//---------------------------------------------------------------------------//
// end DTK_EntityLocalMap.hpp
//---------------------------------------------------------------------------//
