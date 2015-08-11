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
 * \file DTK_ClassicMesh_impl.hpp
 * \author Stuart R. Slattery
 * \brief Mesh manager declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CLASSICMESH_IMPL_HPP
#define DTK_CLASSICMESH_IMPL_HPP

#include "DTK_DBC.hpp"
#include "DTK_MeshTools.hpp"

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class Mesh>
ClassicMesh<Mesh>::ClassicMesh(
    const Teuchos::RCP<MeshManager<Mesh> >& mesh_manager )
    : d_mesh_manager( mesh_manager )
{
    // Allocate arrays.
    int num_blocks = 0;
    if ( Teuchos::nonnull(d_mesh_manager) )
    {
	num_blocks = d_mesh_manager->getNumBlocks();
	d_element_gids.resize( num_blocks );
	d_block_topo.resize( num_blocks );
	d_vertex_g2l.resize( num_blocks );
	d_element_g2l.resize( num_blocks );
	d_vertex_gids.resize( num_blocks );
	d_vertex_coords.resize( num_blocks );
	d_element_conn.resize( num_blocks );
	d_permutation.resize( num_blocks );
    }
    
    // Create views.
    Teuchos::RCP<Mesh> block;
    for ( int b = 0; b < num_blocks; ++b )
    {
	block = d_mesh_manager->getBlock( b );
	d_vertex_gids[b] = MeshTools<Mesh>::verticesView( *block );
	d_vertex_coords[b] = MeshTools<Mesh>::coordsView( *block );
	d_element_gids[b] = MeshTools<Mesh>::elementsView( *block );
	d_element_conn[b] = MeshTools<Mesh>::connectivityView( *block );
	d_permutation[b] = MeshTools<Mesh>::permutationView( *block );
	d_block_topo[b] = createBlockTopology( b );
    }

    // Create global-to-local maps.
    int lid = 0;
    for ( int b = 0; b < num_blocks; ++b )
    {
	lid = 0;
	for ( auto vertex_gid : d_vertex_gids[b] )
	{
	    d_vertex_g2l[b].emplace( vertex_gid, lid );
	    ++lid;
	}
	lid = 0;
	for ( auto element_gid : d_element_gids[b] )
	{
	    d_element_block_map.emplace( element_gid, b );
	    d_element_g2l[b].emplace( element_gid, lid );
	    ++lid;
	}
    }
}

//---------------------------------------------------------------------------//
// Given a block id get its shards topology.
template<class Mesh>
shards::CellTopology
ClassicMesh<Mesh>::createBlockTopology( const int block_id ) const
{
    DTK_REQUIRE( block_id < d_mesh_manager->getNumBlocks() );
    DTK_ElementTopology topo =
	MeshTraits<Mesh>::elementTopology( *getBlock(block_id) );
    int num_nodes = d_permutation[block_id].size();
    const CellTopologyData* topo_data = nullptr;
    switch ( topo )
    {
	case DTK_VERTEX:
	    topo_data = shards::getCellTopologyData<shards::Node>();
	    break;
	case DTK_LINE_SEGMENT:
	    if ( 2 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Line<2> >();
	    else if ( 3 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Line<3> >();
	    break;
	case DTK_TRIANGLE:
	    if ( 3 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Triangle<3> >();
	    else if ( 4 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Triangle<4> >();
	    else if ( 6 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Triangle<6> >();
	    break;
	case DTK_QUADRILATERAL:
	    if ( 4 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Quadrilateral<4> >();
	    else if ( 8 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Quadrilateral<8> >();
	    else if ( 9 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Quadrilateral<9> >();
	    break;
	case DTK_TETRAHEDRON:
	    if ( 4 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Tetrahedron<4> >();
	    else if ( 8 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Tetrahedron<8> >();
	    else if ( 10 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Tetrahedron<10> >();
	    else if ( 11 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Tetrahedron<11> >();
	    break;
	case DTK_PYRAMID:
	    if ( 5 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Pyramid<5> >();
	    else if ( 13 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Pyramid<13> >();
	    else if ( 14 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Pyramid<14> >();
	    break;
	case DTK_WEDGE:
	    if ( 6 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Wedge<6> >();
	    else if ( 15 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Wedge<15> >();
	    else if ( 18 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Wedge<18> >();
	    break;
	case DTK_HEXAHEDRON:
	    if ( 8 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Hexahedron<8> >();
	    else if ( 20 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Hexahedron<20> >();
	    else if ( 27 == num_nodes )
		topo_data = shards::getCellTopologyData<shards::Hexahedron<27> >();
	    break;
    }
    return shards::CellTopology( topo_data );
}

//---------------------------------------------------------------------------//
// Given an element id get its block id.
template<class Mesh>
int ClassicMesh<Mesh>::elementBlockId( const GlobalOrdinal gid ) const
{
    DTK_REQUIRE( d_element_block_map.count(gid) );
    return d_element_block_map.find(gid)->second;
}

//---------------------------------------------------------------------------//
// Given a vertex global id and block id, get the local id in the block.
template<class Mesh>
int ClassicMesh<Mesh>::vertexLocalId( const GlobalOrdinal gid,
				      const int block_id ) const
{
    DTK_REQUIRE( block_id < d_mesh_manager->getNumBlocks() );
    DTK_REQUIRE( d_vertex_g2l[block_id].count(gid) );
    return d_vertex_g2l[block_id].find(gid)->second;
}

//---------------------------------------------------------------------------//
// Given a element global id and block id, get the local id in the block.
template<class Mesh>
int ClassicMesh<Mesh>::elementLocalId( const GlobalOrdinal gid,
				       const int block_id ) const
{
    DTK_REQUIRE( block_id < d_mesh_manager->getNumBlocks() );
    DTK_REQUIRE( d_element_g2l[block_id].count(gid) );
    return d_element_g2l[block_id].find(gid)->second;
}

//---------------------------------------------------------------------------//
// Given an element global id and its block id get the coordinates of the
// element nodes.
template<class Mesh>
Intrepid::FieldContainer<double>
ClassicMesh<Mesh>::getElementNodeCoordinates( const GlobalOrdinal gid,
					      const int block_id ) const
{
    DTK_REQUIRE( block_id < d_mesh_manager->getNumBlocks() );
	
    // Get the element local id.
    int element_lid = elementLocalId( gid, block_id );

    // Allocate the coordinate container.
    int space_dim = dim();
    int num_node = d_permutation[block_id].size();
    Intrepid::FieldContainer<double> coords( 1, num_node, space_dim );

    // Extract the coordinates.
    GlobalOrdinal element_size = d_element_gids[block_id].size();
    GlobalOrdinal vertex_size = d_vertex_gids[block_id].size();
    GlobalOrdinal conn_stride = 0;
    GlobalOrdinal coord_stride = 0;
    GlobalOrdinal vertex_gid = 0;
    int vertex_lid = 0;
    for ( int n = 0; n < num_node; ++n )
    {
	conn_stride = d_permutation[block_id][n] * element_size;
	vertex_gid = d_element_conn[block_id][conn_stride + element_lid];
	vertex_lid = vertexLocalId( vertex_gid, block_id );
	for ( int d = 0; d < space_dim; ++d )
	{
	    coord_stride = d*vertex_size;
	    coords(0,n,d) = d_vertex_coords[block_id][coord_stride + vertex_lid];
	}
    }
    return coords;
}

//---------------------------------------------------------------------------//
// Get the connectivity of an element.
template<class Mesh>
Teuchos::Array<SupportId>
ClassicMesh<Mesh>::getElementConnectivity(
    const GlobalOrdinal gid, const int block_id ) const
{
    DTK_REQUIRE( block_id < d_mesh_manager->getNumBlocks() );
    int element_lid = elementLocalId( gid, block_id );
    int num_node = d_permutation[block_id].size();
    GlobalOrdinal element_size = d_element_gids[block_id].size();
    GlobalOrdinal conn_stride = 0;
    Teuchos::Array<SupportId> conn( num_node );
    for ( int n = 0; n < num_node; ++n )
    {
	conn_stride = d_permutation[block_id][n] * element_size;
	conn[n] = d_element_conn[block_id][conn_stride + element_lid];
    }
    return conn;
}

//---------------------------------------------------------------------------//
// Get the coordinates of a single node.
template<class Mesh>
Teuchos::Array<double>
ClassicMesh<Mesh>::getNodeCoordinates(
    const GlobalOrdinal gid, const int block_id ) const
{
    DTK_REQUIRE( block_id < d_mesh_manager->getNumBlocks() );
    int space_dim = dim();
    int vertex_lid = vertexLocalId( gid, block_id );
    GlobalOrdinal vertex_size = d_vertex_gids[block_id].size();
    Teuchos::Array<double> coords( space_dim );
    for ( int d = 0; d < space_dim; ++d )
    {
	coords[d] = d_vertex_coords[block_id][d*vertex_size + vertex_lid];
    }
    return coords;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_CLASSICMESH_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_ClassicMesh_impl.hpp
//---------------------------------------------------------------------------//

