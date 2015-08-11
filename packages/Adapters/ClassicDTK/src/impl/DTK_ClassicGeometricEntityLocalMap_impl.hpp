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
 * \brief DTK_ClassicGeometricEntityLocalMap_impl.hpp
 * \author Stuart R. Slattery
 * \brief Forward and reverse local mappings for entities.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CLASSICGEOMETRICENTITYLOCALMAP_IMPL_HPP
#define DTK_CLASSICGEOMETRICENTITYLOCALMAP_IMPL_HPP

#include "DTK_ClassicGeometricEntityExtraData.hpp"
#include "DTK_DBC.hpp"

#include "DTK_GeometryTraits.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class Geometry>
ClassicGeometricEntityLocalMap<Geometry>::ClassicGeometricEntityLocalMap()
    : d_inclusion_tol( 1.0e-6 )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Set parameters for mapping.
template<class Geometry>
void ClassicGeometricEntityLocalMap<Geometry>::setParameters(
    const Teuchos::ParameterList& parameters )
{
    if ( parameters.isParameter("Point Inclusion Tolerance") )
    {	    
	d_inclusion_tol = parameters.get<double>("Point Inclusion Tolerance");
    }
}

//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric dimension (volume
// for a 3D entity, area for 2D, and length for 1D).
template<class Geometry>
double ClassicGeometricEntityLocalMap<Geometry>::measure( const Entity& entity ) const
{
    Teuchos::Ptr<Geometry> geometry =
	Teuchos::rcp_dynamic_cast<ClassicGeometricEntityExtraData<Geometry> >(
	    entity.extraData())->d_geometry;
    return GeometryTraits<Geometry>::measure( *geometry );
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
template<class Geometry>
void ClassicGeometricEntityLocalMap<Geometry>::centroid( 
    const Entity& entity, const Teuchos::ArrayView<double>& centroid ) const
{ 
    Teuchos::Ptr<Geometry> geometry =
	Teuchos::rcp_dynamic_cast<ClassicGeometricEntityExtraData<Geometry> >(
	    entity.extraData())->d_geometry;
    centroid.assign(
	GeometryTraits<Geometry>::centroid( *geometry )()
	);
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the parameterized
// point.
template<class Geometry>
bool ClassicGeometricEntityLocalMap<Geometry>::mapToReferenceFrame( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& physical_point,
    const Teuchos::ArrayView<double>& reference_point ) const
{
    reference_point.assign( physical_point );
    return true;
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of an entity.
template<class Geometry>
bool ClassicGeometricEntityLocalMap<Geometry>::checkPointInclusion( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    Teuchos::Ptr<Geometry> geometry =
	Teuchos::rcp_dynamic_cast<ClassicGeometricEntityExtraData<Geometry> >(
	    entity.extraData())->d_geometry;
    Teuchos::Array<double> coords( reference_point );
    return GeometryTraits<Geometry>::pointInGeometry(
	*geometry, coords, d_inclusion_tol );
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
template<class Geometry>
void ClassicGeometricEntityLocalMap<Geometry>::mapToPhysicalFrame( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& physical_point ) const
{
    physical_point.assign( reference_point );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_CLASSICGEOMETRICENTITYLOCALMAP_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_ClassicGeometricEntityLocalMap_impl.hpp
//---------------------------------------------------------------------------//
