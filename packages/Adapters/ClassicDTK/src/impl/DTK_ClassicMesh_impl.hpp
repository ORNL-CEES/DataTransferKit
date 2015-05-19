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
#include "DTK_Classic_MeshTools.hpp"

#include <Shards_BasicTopologies.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template<class Mesh>
ClassicMesh<Mesh>::ClassicMesh(
    const Teuchos::RCP<Classic::MeshManager<Mesh> >& mesh_manager )
    : d_mesh_manager( mesh_manager )
    , d_block_topo( d_mesh_manager->getNumBlocks() )
    , d_vertex_g2l( d_mesh_manager->getNumBlocks() )
    , d_element_g2l( d_mesh_manager->getNumBlocks() )
    , d_vertex_gids( d_mesh_manager->getNumBlocks() )
    , d_vertex_coords( d_mesh_manager->getNumBlocks() )
    , d_element_gids( d_mesh_manager->getNumBlocks() )
    , d_element_conn( d_mesh_manager->getNumBlocks() )
    , d_permutation( d_mesh_manager->getNumBlocks() )
{
    // Create views.
    Teuchos::RCP<Mesh> block;
    for ( int b = 0; b < d_mesh_manager->getNumBlocks(); ++b )
    {
	block = d_mesh_manager->getBlock( b );
	d_block_topo = createBlockTopology( b );
	d_vertex_gids[b] = MeshTools<Mesh>::verticesView( *block );
	d_vertex_coords[b] = MeshTools<Mesh>::coordsView( *block );
	d_element_gids[b] = MeshTools<Mesh>::elementsView( *block );
	d_element_conn[b] = MeshTools<Mesh>::connectivityView( *block );
	d_permutation[b] = MeshTools<Mesh>::permutationView( *block );
    }

    // Create global-to-local maps.
    int lid = 0;
    for ( int b = 0; b < d_mesh_manager->getNumBlocks(); ++b )
    {
	lid = 0;
	for ( auto vertex_gid : d_vertex_gids[b] )
	{
	    d_vertex_g2l[b].emplace( vertex_gid, lid );
	    ++lid
	}
	lid = 0;
	for ( auto element_gid : d_element_gids[b] )
	{
	    d_element_block_map.emplace( element_gid, b );
	    d_element_g2l[b].emplace( element_gid, lid );
	    ++lid
	}
    }
}

//---------------------------------------------------------------------------//
// Given a block id get its shards topology.
shards::CellTopology getBlockTopology( const int block_id ) const
{
    Classic::DTK_Classic_ElementTopology topo =
	Classic::MeshTraits<Mesh>::elementTopology( *getBlock(block_id) );
    int num_nodes = d_permutation[block_id].size();
    const shards::CellTopologyData* topo_data = nullptr;
    switch ( topo )
    {
	case Classic::DTK_Classic_VERTEX:
	    topo_data = getCellTopologyData<shards::Node>();
	    break;
	case Classic::DTK_Classic_LINE_SEGMENT:
	    if ( 2 == num_nodes )
		topo_data = getCellTopologyData<shards::Line<2> >();
	    else if ( 3 == num_nodes )
		topo_data = getCellTopologyData<shards::Line<3> >();
	    break;
	case Classic::DTK_Classic_TRIANGLE:
	    if ( 3 == num_nodes )
		topo_data = getCellTopologyData<shards::Triangle<3> >();
	    else if ( 4 == num_nodes )
		topo_data = getCellTopologyData<shards::Triangle<4> >();
	    else if ( 6 == num_nodes )
		topo_data = getCellTopologyData<shards::Triangle<6> >();
	    break;
	case Classic::DTK_Classic_QUADRILATERAL:
	    if ( 4 == num_nodes )
		topo_data = getCellTopologyData<shards::Quadrilateral<4> >();
	    else if ( 8 == num_nodes )
		topo_data = getCellTopologyData<shards::Quadrilateral<8> >();
	    else if ( 9 == num_nodes )
		topo_data = getCellTopologyData<shards::Quadrilateral<9> >();
	    break;
	case Classic::DTK_Classic_TETRAHEDRON:
	    if ( 4 == num_nodes )
		topo_data = getCellTopologyData<shards::Tetrahedron<4> >();
	    else if ( 8 == num_nodes )
		topo_data = getCellTopologyData<shards::Tetrahedron<8> >();
	    else if ( 10 == num_nodes )
		topo_data = getCellTopologyData<shards::Tetrahedron<10> >();
	    else if ( 11 == num_nodes )
		topo_data = getCellTopologyData<shards::Tetrahedron<11> >();
	    break;
	case Classic::DTK_Classic_PYRAMID:
	    if ( 5 == num_nodes )
		topo_data = getCellTopologyData<shards::Pyramid<5> >();
	    else if ( 13 == num_nodes )
		topo_data = getCellTopologyData<shards::Pyramid<13> >();
	    else if ( 14 == num_nodes )
		topo_data = getCellTopologyData<shards::Pyramid<14> >();
	    break;
	case Classic::DTK_Classic_WEDGE:
	    if ( 6 == num_nodes )
		topo_data = getCellTopologyData<shards::Wedge<6> >();
	    else if ( 15 == num_nodes )
		topo_data = getCellTopologyData<shards::Wedge<15> >();
	    else if ( 18 == num_nodes )
		topo_data = getCellTopologyData<shards::Wedge<18> >();
	    break;
	case Classic::DTK_Classic_HEXAHEDRON:
	    if ( 8 == num_nodes )
		topo_data = getCellTopologyData<shards::Hexahedron<8> >();
	    else if ( 20 == num_nodes )
		topo_data = getCellTopologyData<shards::Hexahedron<20> >();
	    else if ( 27 == num_nodes )
		topo_data = getCellTopologyData<shards::Hexahedron<27> >();
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
    DTK_REQUIRE( d_vertex_l2g[block_id].count(gid) );
    return d_vertex_l2g[block_id].find(gid)->second;
}

//---------------------------------------------------------------------------//
// Given a element global id and block id, get the local id in the block.
template<class Mesh>
int ClassicMesh<Mesh>::elementLocalId( const GlobalOrdinal gid,
				       const int block_id ) const
{
    DTK_REQUIRE( d_element_l2g[block_id].count(gid) );
    return d_element_l2g[block_id].find(gid)->second;
}

//---------------------------------------------------------------------------//
// Given an element global id and its block id get the coordinates of the
// element nodes.
template<class Mesh>
Intrepid::FieldContainer<double>
ClassicMesh<Mesh>::getElementNodeCoordinates( const GlobalOrdinal gid,
					      const int block_id ) const
{
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
	conn_stride = d_permutation[n] * element_size;
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
    int element_lid = elementLocalId( gid, block_id );
    int num_node = d_permutation[block_id].size();
    GlobalOrdinal conn_stride = 0;
    Teuchos::Array<GlobalOrdinal> conn( num_node );
    for ( int n = 0; n < num_node; ++n )
    {
	conn_stride = d_permutation[n] * element_size;
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
    int space_dim = dim();
    int vertex_lid = vertexLocalId( gid, block_id );
    Teuchos::Array<double> coords( space_dim );
    GlobalOrdinal vertex_size = d_vertex_gids[block_id].size();
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

