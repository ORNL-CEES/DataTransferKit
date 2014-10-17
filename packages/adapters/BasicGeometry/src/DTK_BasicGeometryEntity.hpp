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
 * \file DTK_BasicGeometryEntity.hpp
 * \author Stuart R. Slattery
 * \brief BasicGeometryEntity declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_BASICGEOMETRYENTITY_HPP
#define DTK_BASICGEOMETRYENTITY_HPP

#include <iostream>

#include "DTK_Entity.hpp"
#include "DTK_BasicGeometryEntityImpl.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
  \class BasicGeometryEntity
  \brief BasicGeometryEntity interface.
  
  BasicGeometryEntity gives an interface for simple geometries. These objects
  effectivelty define their own EntityLocalMap interface as these functions
  are typically statisfied with analytic expressions for basic geometric
  objects.
*/
//---------------------------------------------------------------------------//
class BasicGeometryEntity : public Entity
{
  public:

    // Default constructor.
    BasicGeometryEntity();

    // Destructor.
    virtual ~BasicGeometryEntity();

    //@{
    //! BasicGeometryEntity interface.
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

#endif // end DTK_BASICGEOMETRYENTITY_HPP

//---------------------------------------------------------------------------//
// end DTK_BasicGeometryEntity.hpp
//---------------------------------------------------------------------------//

