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
 * \brief DTK_EntityFunctionSpace.cpp
 * \author Stuart R. Slattery
 * \brief Discretization functionspace interface.
 */
//---------------------------------------------------------------------------//

#include "DTK_EntityFunctionSpace.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
EntityFunctionSpace::EntityFunctionSpace()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Destructor.
EntityFunctionSpace::~EntityFunctionSpace()
{ /* ... */ }

//---------------------------------------------------------------------------//
// Return the entity measure with respect to the parameteric dimension (volume
// for a 3D entity, area for 2D, and length for 1D). 
double EntityFunctionSpace::measure( const Entity& entity ) const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
    return -1.0;
}

//---------------------------------------------------------------------------//
// Return the centroid of the entity.
void EntityFunctionSpace::centroid(
    const Entity& entity,
    Teuchos::ArrayView<const double>& centroid ) const
{ 
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Perform a safeguard check for mapping a point to the reference space
// of an entity using the given tolerance. 
void EntityFunctionSpace::safeguardMapToReferenceFrame(
    const Entity& entity,
    const Teuchos::ParameterList& parameters,
    const Teuchos::ArrayView<const double>& point,
    MappingStatus& status ) const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Map a point to the reference space of an entity. Return the parameterized point.
void EntityFunctionSpace::mapToReferenceFrame( 
    const Entity& entity,
    const Teuchos::ParameterList& parameters,
    const Teuchos::ArrayView<const double>& point,
    const Teuchos::ArrayView<double>& reference_point,
    MappingStatus& status ) const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Determine if a reference point is in the parameterized space of an entity.
bool EntityFunctionSpace::checkPointInclusion( 
    const Entity& entity,
    const Teuchos::ParameterList& parameters,
    const Teuchos::ArrayView<const double>& reference_point ) const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
    return false;
}

//---------------------------------------------------------------------------//
// Map a reference point to the physical space of an entity.
void EntityFunctionSpace::mapToPhysicalFrame( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    const Teuchos::ArrayView<double>& point ) const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//
// Given an entity and a reference point, evaluate the functional support
// of the entity at that point. 
void EntityFunctionSpace::evaluateShapeFunction( 
    const Entity& entity,
    const Teuchos::ArrayView<const double>& reference_point,
    Teuchos::Array<EntityId>& dof_ids,
    Kokkos::View<double>& values ) const
{
    DTK_REMEMBER( bool not_implemented = true );
    DTK_INSIST( !not_implemented );
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_EntityFunctionSpace.cpp
//---------------------------------------------------------------------------//
