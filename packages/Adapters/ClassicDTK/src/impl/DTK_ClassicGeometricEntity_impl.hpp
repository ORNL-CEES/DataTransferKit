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
 * \file DTK_ClassicGeometricEntity_impl.hpp
 * \author Stuart R. Slattery
 * \brief ClassicGeometricEntity definition
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CLASSIC_GEOMETRICENTITY_IMPL_HPP
#define DTK_CLASSIC_GEOMETRICENTITY_IMPL_HPP

#include <limits>

#include "DTK_DBC.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Default constructor.
template<class Geometry>
ClassicGeometricEntity<Geometry>::ClassicGeometricEntity(
    const Teuchos::Ptr<Geometry>& geometry,
    const EntityId global_id,
    const int owner_rank )
{
    this->b_entity_impl = Teuchos::rcp(
	new ClassicGeometricEntityImpl(geometry,global_id,rank) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the measure of the entity.
 *
 * \return Return the measure of the entity.
 */
template<class Geometry>
double ClassicGeometricEntity<Geometry>::measure() const
{
    DTK_REQUIRE( Teuchos::nonnull(b_entity_impl) );
    return Teuchos::rcp_dynamic_cast<ClassicGeometricEntityImpl>(
	this->b_entity_impl)->measure();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the centroid of the entity.
 *
 * \return The centroid coordinates.
 */
template<class Geometry>
void ClassicGeometricEntity<Geometry>::centroid( 
    const Teuchos::ArrayView<double>& centroid ) const
{
    DTK_REQUIRE( Teuchos::nonnull(b_entity_impl) );
    Teuchos::rcp_dynamic_cast<ClassicGeometricEntityImpl>(
	this->b_entity_impl)->centroid(centroid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map a point to the reference space of an entity. Return the
 */
template<class Geometry>
bool ClassicGeometricEntity<Geometry>::mapToReferenceFrame( 
    const Teuchos::ArrayView<const double>& point,
    const Teuchos::ArrayView<double>& reference_point ) const
{
    DTK_REQUIRE( Teuchos::nonnull(b_entity_impl) );
    return Teuchos::rcp_dynamic_cast<ClassicGeometricEntityImpl>(
	this->b_entity_impl)->mapToReferenceFrame(point,reference_point);
}

//---------------------------------------------------------------------------//
/*!  
 * \brief Determine if a reference point is in the parameterized space of
 * an entity.
 */
template<class Geometry>
bool ClassicGeometricEntity<Geometry>::checkPointInclusion( 
    const double tolerance,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    DTK_REQUIRE( Teuchos::nonnull(b_entity_impl) );
    return Teuchos::rcp_dynamic_cast<ClassicGeometricEntityImpl>(
	this->b_entity_impl)->checkPointInclusion(tolerance,reference_point);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map a reference point to the physical space of an entity.
 */
template<class Geometry>
void ClassicGeometricEntity<Geometry>::mapToPhysicalFrame( 
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& point ) const
{
    DTK_REQUIRE( Teuchos::nonnull(b_entity_impl) );
    Teuchos::rcp_dynamic_cast<ClassicGeometricEntityImpl>(
	this->b_entity_impl)->mapToPhysicalFrame(reference_point,point);
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_CLASSIC_GEOMETRICENTITY_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_ClassicGeometricEntity_impl.hpp
//---------------------------------------------------------------------------//

