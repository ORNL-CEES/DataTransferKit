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
 * \brief DTK_LibmeshEntityImpl.cpp
 * \author Stuart R. Slattery
 * \brief Libmesh entity implementation.
 */
//---------------------------------------------------------------------------//

#include <limits>
#include <vector>

#include "DTK_LibmeshEntityImpl.hpp"
#include <DTK_DBC.hpp>

#include <libmesh/elem.h>
#include <libmesh/node.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// libMesh::Elem specializations.
//---------------------------------------------------------------------------//
// Constructor.
template <>
LibmeshEntityImpl<libMesh::Elem>::LibmeshEntityImpl(
    const Teuchos::Ptr<libMesh::Elem> &libmesh_geom,
    const Teuchos::Ptr<libMesh::MeshBase> &libmesh_mesh,
    const Teuchos::Ptr<LibmeshAdjacencies> &adjacencies )
    : d_extra_data( new LibmeshEntityExtraData<libMesh::Elem>( libmesh_geom ) )
    , d_mesh( libmesh_mesh )
    , d_adjacencies( adjacencies )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
template <>
DataTransferKit::EntityId LibmeshEntityImpl<libMesh::Elem>::id() const
{
    DTK_REQUIRE( d_extra_data->d_libmesh_geom->valid_id() );
    return d_extra_data->d_libmesh_geom->id();
}

//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
template <>
int LibmeshEntityImpl<libMesh::Elem>::ownerRank() const
{
    DTK_REQUIRE( d_extra_data->d_libmesh_geom->valid_processor_id() );
    return d_extra_data->d_libmesh_geom->processor_id();
}

//---------------------------------------------------------------------------//
// Return the topological dimension of the entity.
template <>
int LibmeshEntityImpl<libMesh::Elem>::topologicalDimension() const
{
    return d_extra_data->d_libmesh_geom->dim();
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
template <>
int LibmeshEntityImpl<libMesh::Elem>::physicalDimension() const
{
    return d_mesh->mesh_dimension();
}

//---------------------------------------------------------------------------//
// Get the extra data on the entity.
template <>
Teuchos::RCP<DataTransferKit::EntityExtraData>
LibmeshEntityImpl<libMesh::Elem>::extraData() const
{
    return d_extra_data;
}

//---------------------------------------------------------------------------//
// Return the Cartesian bounding box around an entity.
template <>
void LibmeshEntityImpl<libMesh::Elem>::boundingBox(
    Teuchos::Tuple<double, 6> &bounds ) const
{
    unsigned int num_nodes = d_extra_data->d_libmesh_geom->n_nodes();
    int space_dim = this->physicalDimension();
    double max = std::numeric_limits<double>::max();
    bounds = Teuchos::tuple( max, max, max, -max, -max, -max );
    for ( unsigned int n = 0; n < num_nodes; ++n )
    {
        const libMesh::Point &node = d_extra_data->d_libmesh_geom->point( n );
        for ( int d = 0; d < space_dim; ++d )
        {
            bounds[d] = std::min( bounds[d], node( d ) );
            bounds[d + 3] = std::max( bounds[d + 3], node( d ) );
        }
    }
    for ( int d = space_dim; d < 3; ++d )
    {
        bounds[d] = -max;
        bounds[d + 3] = max;
    }
}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
template <>
bool LibmeshEntityImpl<libMesh::Elem>::inBlock( const int block_id ) const
{
    return ( block_id == d_extra_data->d_libmesh_geom->subdomain_id() );
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
template <>
bool LibmeshEntityImpl<libMesh::Elem>::onBoundary( const int boundary_id ) const
{
    bool on_boundary = false;
    int n_sides = d_extra_data->d_libmesh_geom->n_sides();
    for ( int s = 0; s < n_sides; ++s )
    {
        on_boundary = d_mesh->get_boundary_info().has_boundary_id(
            d_extra_data->d_libmesh_geom.getRawPtr(), s, boundary_id );
        if ( on_boundary )
        {
            break;
        }
    }
    return on_boundary;
}

//---------------------------------------------------------------------------//
// Provide a verbose description of the object.
template <>
void LibmeshEntityImpl<libMesh::Elem>::describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel /*verb_level*/ ) const
{
    out << std::endl;
    out << "---" << std::endl;
    out << "LibMesh Element" << std::endl;
    out << "Owner rank: " << ownerRank() << std::endl;
    out << d_extra_data->d_libmesh_geom->get_info() << std::endl;
    out << "---" << std::endl;
}

//---------------------------------------------------------------------------//
// libMesh::Node specializations.
//---------------------------------------------------------------------------//
// Constructor.
template <>
LibmeshEntityImpl<libMesh::Node>::LibmeshEntityImpl(
    const Teuchos::Ptr<libMesh::Node> &libmesh_geom,
    const Teuchos::Ptr<libMesh::MeshBase> &libmesh_mesh,
    const Teuchos::Ptr<LibmeshAdjacencies> &adjacencies )
    : d_extra_data( new LibmeshEntityExtraData<libMesh::Node>( libmesh_geom ) )
    , d_mesh( libmesh_mesh )
    , d_adjacencies( adjacencies )
{ /* ... */
}

//---------------------------------------------------------------------------//
// Get the unique global identifier for the entity.
template <>
DataTransferKit::EntityId LibmeshEntityImpl<libMesh::Node>::id() const
{
    DTK_REQUIRE( d_extra_data->d_libmesh_geom->valid_id() );
    return d_extra_data->d_libmesh_geom->id();
}

//---------------------------------------------------------------------------//
// Get the parallel rank that owns the entity.
template <>
int LibmeshEntityImpl<libMesh::Node>::ownerRank() const
{
    DTK_REQUIRE( d_extra_data->d_libmesh_geom->valid_processor_id() );
    return d_extra_data->d_libmesh_geom->processor_id();
}

//---------------------------------------------------------------------------//
// Return the topological dimension of the entity.
template <>
int LibmeshEntityImpl<libMesh::Node>::topologicalDimension() const
{
    return 0;
}

//---------------------------------------------------------------------------//
// Return the physical dimension of the entity.
template <>
int LibmeshEntityImpl<libMesh::Node>::physicalDimension() const
{
    return d_mesh->mesh_dimension();
}

//---------------------------------------------------------------------------//
// Get the extra data on the entity.
template <>
Teuchos::RCP<DataTransferKit::EntityExtraData>
LibmeshEntityImpl<libMesh::Node>::extraData() const
{
    return d_extra_data;
}

//---------------------------------------------------------------------------//
// Return the Cartesian bounding box around an entity.
template <>
void LibmeshEntityImpl<libMesh::Node>::boundingBox(
    Teuchos::Tuple<double, 6> &bounds ) const
{
    int space_dim = this->physicalDimension();
    for ( int d = 0; d < space_dim; ++d )
    {
        bounds[d] = ( *( d_extra_data->d_libmesh_geom ) )( d );
        bounds[d + 3] = ( *( d_extra_data->d_libmesh_geom ) )( d );
    }
}

//---------------------------------------------------------------------------//
// Determine if an entity is in the block with the given id.
template <>
bool LibmeshEntityImpl<libMesh::Node>::inBlock( const int block_id ) const
{
    Teuchos::Array<Teuchos::Ptr<libMesh::Elem>> node_elems;
    d_adjacencies->getLibmeshAdjacencies( d_extra_data->d_libmesh_geom,
                                          node_elems );
    for ( auto &elem : node_elems )
    {
        if ( block_id == elem->subdomain_id() )
        {
            return true;
        }
    }

    return false;
}

//---------------------------------------------------------------------------//
// Determine if an entity is on the boundary with the given id.
template <>
bool LibmeshEntityImpl<libMesh::Node>::onBoundary( const int boundary_id ) const
{
    return d_mesh->get_boundary_info().has_boundary_id(
        d_extra_data->d_libmesh_geom.getRawPtr(), boundary_id );
}

//---------------------------------------------------------------------------//
// Provide a verbose description of the object.
template <>
void LibmeshEntityImpl<libMesh::Node>::describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel /*verb_level*/ ) const
{
    out << std::endl;
    out << "---" << std::endl;
    out << "LibMesh Node" << std::endl;
    out << "Id: " << id() << std::endl;
    out << "Owner rank: " << ownerRank() << std::endl;
    out << d_extra_data->d_libmesh_geom->get_info() << std::endl;
    out << "---" << std::endl;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_LibmeshEntityImpl.cpp
//---------------------------------------------------------------------------//
