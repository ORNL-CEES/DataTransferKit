//---------------------------------------------------------------------------//
/*
  Copyright (c) 2014, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the Oak Ridge National Laboratory nor the
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
 * \brief DTK_IntrepidCellLocalMap.hpp
 * \author Stuart R. Slattery
 * \brief Forward and reverse local mappings for IntrepidCells.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTREPIDCELLLOCALMAP_HPP
#define DTK_INTREPIDCELLLOCALMAP_HPP

#include <unordered_map>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ArrayView.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_FieldContainer.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class IntrepidcellLocalMap
  \brief A stateless class of IntrepidCell wrappers for implementing
  EntityLocalMap for element-level entities.
*/
//---------------------------------------------------------------------------//
class IntrepidCellLocalMap
{
  public:

    /*!
     * \brief Constructor.
     */
    IntrepidCellLocalMap() { /* ... */ }

    /*!
     * \brief Destructor.
     */
    ~IntrepidCellLocalMap() { /* ... */ }

    /*!
     * \brief Return the entity measure with respect to the parameteric
     * dimension (volume for a 3D entity, area for 2D, and length for 1D).
     * \param entity Compute the measure for this entity.
     * \return The measure of the entity.
     */
    static double measure( const shards::CellTopology& entity_topo,
			   const Intrepid::FieldContainer<double>& entity_coords );

    /*!
     * \brief Return the centroid of the entity.
     * \param centroid A view of the centroid coordinates. This view will
     * be allocated. Assign a view of your centroid to this view.
     */
    static void centroid( const shards::CellTopology& entity_topo,
			  const Intrepid::FieldContainer<double>& entity_coords,
			  const Teuchos::ArrayView<double>& centroid );

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
    static bool mapToReferenceFrame( 
	const shards::CellTopology& entity_topo,
	const Intrepid::FieldContainer<double>& entity_coords,
	const Teuchos::ArrayView<const double>& point,
	const Teuchos::ArrayView<double>& reference_point );

    /*!  
     * \brief Determine if a reference point is in the parameterized space of
     * an entity.
     * \param entity Perfrom the mapping for this entity.
     * \param parameters Parameters to be used for the point inclusion check.
     * \param reference_point A view into an array of size physicalDimension()
     * containing the reference coordinates of the mapped point.
     * \return True if the point is in the reference space, false if not.
     */
    static bool checkPointInclusion( 
	const shards::CellTopology& entity_topo,
	const Teuchos::ArrayView<const double>& reference_point,
	const double tolerance );

    /*!
     * \brief (Forward Map) Map a reference point to the physical space of an
     * entity. 
     * \param entity Perfrom the mapping for this entity.
     * \param reference_point A view into an array of size physicalDimension()
     * containing the reference coordinates of the mapped point.
     * \param A view into an array of size physicalDimension() to write
     * the coordinates of physical point.
     */
    static void mapToPhysicalFrame( 
	const shards::CellTopology& entity_topo,
	const Intrepid::FieldContainer<double>& entity_coords,
	const Teuchos::ArrayView<const double>& reference_point,
	const Teuchos::ArrayView<double>& point );
};

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_INTREPIDCELLLOCALMAP_HPP

//---------------------------------------------------------------------------//
// end DTK_IntrepidCellLocalMap.hpp
//---------------------------------------------------------------------------//
