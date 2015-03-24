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
 * \brief DTK_MoabEntityLocalMap.hpp
 * \author Stuart R. Slattery
 * \brief Forward and reverse local mappings for entities.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOABENTITYLOCALMAP_HPP
#define DTK_MOABENTITYLOCALMAP_HPP

#include "DTK_EntityLocalMap.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ArrayView.hpp>

#include <MBParallelComm.hpp>
#include <moab/ElemEvaluator.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class MoabEntityLocalMap
  \brief Moab mesh forward and reverse local map implementation.
*/
//---------------------------------------------------------------------------//
class MoabEntityLocalMap : public EntityLocalMap
{
  public:

    /*!
     * \brief Constructor.
     */
    MoabEntityLocalMap( const Teuchos::RCP<moab::ParallelComm>& moab_mesh );

    /*!
     * \brief Return the entity measure with respect to the parameteric
     * dimension (volume for a 3D entity, area for 2D, and length for 1D).
     * \param entity Compute the measure for this entity.
     * \return The measure of the entity.
     */
    double measure( const Entity& entity ) const override;

    /*!
     * \brief Return the centroid of the entity.
     * \param centroid A view of the centroid coordinates. This view will
     * be allocated. Assign a view of your centroid to this view.
     */
    void centroid( const Entity& entity,
		   const Teuchos::ArrayView<double>& centroid ) const override;

    /*!
     * \brief (Safeguard the reverse map) Perform a safeguard check for
     * mapping a point to the reference space of an entity using the given
     * tolerance.
     * \param entity Perfrom the mapping for this entity.
     * \param parameters Parameters to be used for the safeguard check.
     * \param point A view into an array of size physicalDimension() containing
     * the coordinates of the point to map.
     * \param status A status object indicating the results of the safeguard
     * check.
     * \return Return true if it is safe to map to the reference frame.
     */
    bool isSafeToMapToReferenceFrame(
	const Entity& entity,
	const Teuchos::ArrayView<const double>& point,
	const Teuchos::RCP<MappingStatus>& status = Teuchos::null ) const override;

    /*!
     * \brief (Reverse Map) Map a point to the reference space of an
     * entity. Return the parameterized point.
     * \param entity Perfrom the mapping for this entity.
     * \param parameters Parameters to be used for the mapping procedure.
     * \param  A view into an array of size physicalDimension() containing
     * the coordinates of the point to map.
     * \param reference_point A view into an array of size physicalDimension()
     * to write the reference coordinates of the mapped point.
     * \param status A status object indicating the results of the mapping
     * procedure.
     * \return Return true if the map to reference frame succeeded.
     */
    bool mapToReferenceFrame( 
	const Entity& entity,
	const Teuchos::ArrayView<const double>& point,
	const Teuchos::ArrayView<double>& reference_point,
	const Teuchos::RCP<MappingStatus>& status = Teuchos::null ) const override;

    /*!  
     * \brief Determine if a reference point is in the parameterized space of
     * an entity.
     * \param entity Perfrom the mapping for this entity.
     * \param parameters Parameters to be used for the point inclusion check.
     * \param reference_point A view into an array of size physicalDimension()
     * containing the reference coordinates of the mapped point.
     * \return True if the point is in the reference space, false if not.
     */
    bool checkPointInclusion( 
	const Entity& entity,
	const Teuchos::ArrayView<const double>& reference_point ) const override;

    /*!
     * \brief (Forward Map) Map a reference point to the physical space of an
     * entity. 
     * \param entity Perfrom the mapping for this entity.
     * \param reference_point A view into an array of size physicalDimension()
     * containing the reference coordinates of the mapped point.
     * \param A view into an array of size physicalDimension() to write
     * the coordinates of physical point.
     */
    void mapToPhysicalFrame( 
	const Entity& entity,
	const Teuchos::ArrayView<const double>& reference_point,
	const Teuchos::ArrayView<double>& point ) const override;

    /*!
     * \brief Compute the normal on a face (3D) or edge (2D) at a given
     * reference point. A default implementation is provided using a finite
     * difference scheme.
     * \param entity Compute the normal for this entity.
     * \param reference_point Compute the normal at this reference point.
     * \param normal A view into an array of size physicalDimension() to write
     * the normal.
     */
    void normalAtReferencePoint( 
        const Entity& entity,
        const Teuchos::ArrayView<double>& reference_point,
        const Teuchos::ArrayView<double>& normal ) const override;

  private:

    // Cache an entity in the evaluator.
    void cacheEntity( const Entity& entity ) const;

    // Get the parameteric center of an entity.
    void parametricCenter( const Entity& entity,
			   Teuchos::Array<double>& center ) const;

  private:

    // Moab mesh.
    Teuchos::RCP<moab::ParallelComm> d_moab_mesh;

    // Moab element evaluator.
    mutable Teuchos::RCP<moab::ElemEvaluator> d_moab_evaluator;
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_MOABENTITYLOCALMAP_HPP

//---------------------------------------------------------------------------//
// end DTK_MoabEntityLocalMap.hpp
//---------------------------------------------------------------------------//
