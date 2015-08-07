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
 * \file DTK_ClassicMeshElementImpl_impl.hpp
 * \author Stuart R. Slattery
 * \brief ClassicMeshElementImpl declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CLASSIC_MESHELEMENTENTITYIMPL_IMPL_HPP
#define DTK_CLASSIC_MESHELEMENTENTITYIMPL_IMPL_HPP

#include "Intrepid_FieldContainer.hpp"

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
template<class Mesh>
ClassicMeshElementImpl<Mesh>::ClassicMeshElementImpl(
    const Teuchos::Ptr<ClassicMesh<Mesh> >& mesh,
    const EntityId global_id,
    const int block_id )
    : d_mesh( mesh )
    , d_id( global_id )
{
    d_extra_data = Teuchos::rcp( new ClassicMeshElementExtraData(block_id) );
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
template<class Mesh>
EntityId ClassicMeshElementImpl<Mesh>::id() const
{
    return d_id;
}
    
//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
template<class Mesh>
int ClassicMeshElementImpl<Mesh>::ownerRank() const
{
    return d_mesh->comm()->getRank();
}

//---------------------------------------------------------------------------//
// Return the topological dimension of the entity.
template<class Mesh>
int ClassicMeshElementImpl<Mesh>::topologicalDimension() const
{
    Classic::DTK_Classic_ElementTopology topo =
	Classic::MeshTraits<Mesh>::elementTopology(
	    *d_mesh->getBlock(d_extra_data->d_block_id) );
    int dim = 0;
    switch ( topo )
    {
	case Classic::DTK_Classic_VERTEX:
	    dim = 0;
	    break;
	case Classic::DTK_Classic_LINE_SEGMENT:
	    dim = 1;
	    break;
	case Classic::DTK_Classic_TRIANGLE:
	    dim = 2;
	    break;
	case Classic::DTK_Classic_QUADRILATERAL:
	    dim = 2;
	    break;
	case Classic::DTK_Classic_TETRAHEDRON:
	    dim = 3;
	    break;
	case Classic::DTK_Classic_PYRAMID:
	    dim = 3;
	    break;
	case Classic::DTK_Classic_WEDGE:
	    dim = 3;
	    break;
	case Classic::DTK_Classic_HEXAHEDRON:
	    dim = 3;
	    break;
	default:
	    dim = -1;
	    break;
    }
    return dim;
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
template<class Mesh>
int ClassicMeshElementImpl<Mesh>::physicalDimension() const
{
    return d_mesh->dim();
}

//---------------------------------------------------------------------------//
// Return the Cartesian bounding box around an entity.
template<class Mesh>
void ClassicMeshElementImpl<Mesh>::boundingBox(
    Teuchos::Tuple<double,6>& bounds ) const
{
    Intrepid::FieldContainer<double> coords =
	d_mesh->getElementNodeCoordinates( d_id, d_extra_data->d_block_id );
    int num_nodes = coords.dimension( 1 );
    int space_dim = coords.dimension( 2 );
    double max = std::numeric_limits<double>::max();
    bounds = Teuchos::tuple( max, max, max, -max, -max, -max );
    for ( int n = 0; n < num_nodes; ++n )
    {
	for ( int d = 0; d < space_dim; ++d )
	{
	    bounds[d] = std::min( bounds[d], coords(0,n,d) );
	    bounds[d+3] = std::max( bounds[d+3], coords(0,n,d) );
	}
    }
}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
template<class Mesh>
bool ClassicMeshElementImpl<Mesh>::inBlock( const int block_id ) const
{
    return false;
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
template<class Mesh>
bool ClassicMeshElementImpl<Mesh>::onBoundary( const int boundary_id ) const
{
    return false;
}

//---------------------------------------------------------------------------//
// Get the extra data on the entity.
template<class Mesh>
Teuchos::RCP<EntityExtraData> ClassicMeshElementImpl<Mesh>::extraData() const
{
    return d_extra_data;
}

//---------------------------------------------------------------------------//
// Provide a one line description of the object.
template<class Mesh>
std::string ClassicMeshElementImpl<Mesh>::description() const
{
    return std::string("DTK Classic Mesh Element");
}

//---------------------------------------------------------------------------//
// Provide a verbose description of the object.
template<class Mesh>
void ClassicMeshElementImpl<Mesh>::describe(
    Teuchos::FancyOStream& out,
    const Teuchos::EVerbosityLevel verb_level ) const
{
    out << description();
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_CLASSIC_MESHELEMENTENTITYIMPL_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_ClassicMeshElementImpl_impl.hpp
//---------------------------------------------------------------------------//

