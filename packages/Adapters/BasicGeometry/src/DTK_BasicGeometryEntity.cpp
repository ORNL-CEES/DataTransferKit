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
 * \file DTK_BasicGeometryEntity.cpp
 * \author Stuart R. Slattery
 * \brief BasicGeometryEntity definition
 */
//---------------------------------------------------------------------------//

#include <limits>

#include "DTK_DBC.hpp"
#include "DTK_BasicGeometryEntity.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Default constructor.
BasicGeometryEntity::BasicGeometryEntity()
{
    this->b_entity_impl = Teuchos::rcp( new BasicGeometryEntityImpl() );
}

//---------------------------------------------------------------------------//
// Destructor.
BasicGeometryEntity::~BasicGeometryEntity()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the measure of the entity.
 *
 * \return Return the measure of the entity.
 */
double BasicGeometryEntity::measure() const
{
    return Teuchos::rcp_dynamic_cast<BasicGeometryEntityImpl>(
	this->b_entity_impl)->measure();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the centroid of the entity.
 *
 * \return The centroid coordinates.
 */
void BasicGeometryEntity::centroid( 
    const Teuchos::ArrayView<double>& centroid ) const
{
    Teuchos::rcp_dynamic_cast<BasicGeometryEntityImpl>(
	this->b_entity_impl)->centroid(centroid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Safeguard the reverse map.
 */
bool BasicGeometryEntity::isSafeToMapToReferenceFrame(
    const Teuchos::ArrayView<const double>& point ) const
{
    return Teuchos::rcp_dynamic_cast<BasicGeometryEntityImpl>(
	this->b_entity_impl)->isSafeToMapToReferenceFrame(point);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map a point to the reference space of an entity. Return the
 */
bool BasicGeometryEntity::mapToReferenceFrame( 
    const Teuchos::ArrayView<const double>& point,
    const Teuchos::ArrayView<double>& reference_point ) const
{
    return Teuchos::rcp_dynamic_cast<BasicGeometryEntityImpl>(
	this->b_entity_impl)->mapToReferenceFrame(point,reference_point);
}

//---------------------------------------------------------------------------//
/*!  
 * \brief Determine if a reference point is in the parameterized space of
 * an entity.
 */
bool BasicGeometryEntity::checkPointInclusion( 
    const double tolerance,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    return Teuchos::rcp_dynamic_cast<BasicGeometryEntityImpl>(
	this->b_entity_impl)->checkPointInclusion(tolerance,reference_point);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map a reference point to the physical space of an entity.
 */
void BasicGeometryEntity::mapToPhysicalFrame( 
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& point ) const
{
    Teuchos::rcp_dynamic_cast<BasicGeometryEntityImpl>(
	this->b_entity_impl)->mapToPhysicalFrame(reference_point,point);
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_BasicGeometryEntity.cpp
//---------------------------------------------------------------------------//

