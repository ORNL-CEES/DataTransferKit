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
 * \file DTK_ClassicGeometricEntityImpl_impl.hpp
 * \author Stuart R. Slattery
 * \brief ClassicGeometricEntityImpl declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CLASSIC_GEOMETRICENTITYENTITYIMPL_IMPL_HPP
#define DTK_CLASSIC_GEOMETRICENTITYENTITYIMPL_IMPL_HPP

#include "DTK_Classic_GeometryTraits.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
template<class Geometry>
ClassicGeometricEntityImpl<Geometry>::ClassicGeometricEntityImpl(
    const Teuchos::Ptr<Geometry>& geometry,
    const EntityId global_id,
    const int owner_rank )
    : d_id( global_id )
    , d_owner_rank( owner_rank )
{
    d_extra_data = Teuchos::rcp(
	new ClassicGeometricEntityExtraData<Geometry>(geometry) );
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
template<class Geometry>
EntityId ClassicGeometricEntityImpl<Geometry>::id() const
{
    return d_id;
}
    
//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
template<class Geometry>
int ClassicGeometricEntityImpl<Geometry>::ownerRank() const
{
    return d_owner_rank;
}

//---------------------------------------------------------------------------//
// Return the topological dimension of the entity.
template<class Geometry>
int ClassicGeometricEntityImpl<Geometry>::topologicalDimension() const
{
    return
	Classic::GeometryTraits<Geometry>::dim( *(d_extra_data->d_geometry) );
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
template<class Geometry>
int ClassicGeometricEntityImpl<Geometry>::physicalDimension() const
{
    return
	Classic::GeometryTraits<Geometry>::dim( *(d_extra_data->d_geometry) );
}

//---------------------------------------------------------------------------//
// Return the Cartesian bounding box around an entity.
template<class Geometry>
void ClassicGeometricEntityImpl<Geometry>::boundingBox(
    Teuchos::Tuple<double,6>& bounds ) const
{
    bounds = Classic::GeometryTraits<Geometry>::boundingBox(
	*(d_extra_data->d_geometry) ).getBounds();
}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
template<class Geometry>
bool ClassicGeometricEntityImpl<Geometry>::inBlock( const int block_id ) const
{
    return false;
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
template<class Geometry>
bool ClassicGeometricEntityImpl<Geometry>::onBoundary( const int boundary_id ) const
{
    return false;
}

//---------------------------------------------------------------------------//
// Get the extra data on the entity.
template<class Geometry>
Teuchos::RCP<EntityExtraData> ClassicGeometricEntityImpl<Geometry>::extraData() const
{
    return d_extra_data;
}

//---------------------------------------------------------------------------//
// Provide a one line description of the object.
template<class Geometry>
std::string ClassicGeometricEntityImpl<Geometry>::description() const
{
    return std::string("DTK Classic Geometric Entity");
}

//---------------------------------------------------------------------------//
// Provide a verbose description of the object.
template<class Geometry>
void ClassicGeometricEntityImpl<Geometry>::describe(
    Teuchos::FancyOStream& out,
    const Teuchos::EVerbosityLevel verb_level ) const
{
    out << description();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_CLASSIC_GEOMETRICENTITYENTITYIMPL_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_ClassicGeometricEntityImpl_impl.hpp
//---------------------------------------------------------------------------//

