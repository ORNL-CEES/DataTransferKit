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
 * \file DTK_RCB_def.hpp
 * \author Stuart R. Slattery
 * \brief Wrapper definition for Zoltan recursive coordinate bisectioning.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "DTK_MeshTools.hpp"
#include "DTK_Assertion.hpp"
#include "DTK_CommIndexer.hpp"
#include "DataTransferKit_config.hpp"

#include <mpi.h>

#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Tuple.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param comm The communicator over which to build the RCB partitioning.
 *
 * \param mesh_manager The mesh to be partitioned with RCB. A null RCP is
 * valid here as the mesh may or may not exist on all of the processes we want
 * to repartition it to.
 *
 * \param dimension The dimension of the RCB space.
 */
template<class Mesh>
RCB<Mesh>::RCB( const RCP_Comm& comm, const RCP_MeshManager& mesh_manager, 
		const int dimension )
    : d_mesh_manager( mesh_manager )
    , d_dimension( dimension )
{
    // Get the raw MPI communicator.
    Teuchos::RCP< const Teuchos::MpiComm<int> > mpi_comm = 
	Teuchos::rcp_dynamic_cast< const Teuchos::MpiComm<int> >( comm );
    Teuchos::RCP< const Teuchos::OpaqueWrapper<MPI_Comm> > opaque_comm = 
	mpi_comm->getRawMpiComm();
    MPI_Comm raw_comm = (*opaque_comm)();

    // Create the Zoltan object.
    d_zz = Zoltan_Create( raw_comm );

    // General parameters.
    Zoltan_Set_Param( d_zz, "DEBUG_LEVEL", "0" );
    Zoltan_Set_Param( d_zz, "LB_METHOD", "RCB" );
    Zoltan_Set_Param( d_zz, "NUM_GID_ENTRIES", "1" ); 
    Zoltan_Set_Param( d_zz, "NUM_LID_ENTRIES", "1" );
    Zoltan_Set_Param( d_zz, "DEBUG_PROCESSOR", "0" );
    Zoltan_Set_Param( d_zz, "OBJ_WEIGHT_DIM", "0" );
    Zoltan_Set_Param( d_zz, "EDGE_WEIGHT_DIM", "0" );
    Zoltan_Set_Param( d_zz, "RETURN_LISTS", "ALL" );

    // RCB parameters.
    Zoltan_Set_Param( d_zz, "RCB_OUTPUT_LEVEL", "0" );
    Zoltan_Set_Param( d_zz, "RCB_RECTILINEAR_BLOCKS", "0" );
    Zoltan_Set_Param( d_zz, "KEEP_CUTS", "1" );
    Zoltan_Set_Param( d_zz, "AVERAGE_CUTS", "1" );
    Zoltan_Set_Param( d_zz, "RCB_LOCK_DIRECTIONS", "1" );
    Zoltan_Set_Param( d_zz, "RCB_SET_DIRECTIONS", "1" );

    // Register static functions.
    Zoltan_Set_Num_Obj_Fn( d_zz, getNumberOfObjects, &d_mesh_manager );
    Zoltan_Set_Obj_List_Fn( d_zz, getObjectList, &d_mesh_manager );
    Zoltan_Set_Num_Geom_Fn( d_zz, getNumGeometry, &d_dimension );
    Zoltan_Set_Geom_Multi_Fn( d_zz, getGeometryList, &d_mesh_manager );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor. Zoltan memory deallocation happens here and only here.
 */
template<class Mesh>
RCB<Mesh>::~RCB()
{
    // Zoltan cleanup.
    Zoltan_LB_Free_Part( &d_import_global_ids, &d_import_local_ids, 
			 &d_import_procs, &d_import_to_part );
    Zoltan_LB_Free_Part( &d_export_global_ids, &d_export_local_ids, 
			 &d_export_procs, &d_export_to_part );
    Zoltan_Destroy( &d_zz );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute RCB partitioning of the mesh.
 */
template<class Mesh>
void RCB<Mesh>::partition()
{
    // Run zoltan partitioning.
    rememberValue( int zoltan_error );
#if HAVE_DTK_DBC
    zoltan_error = Zoltan_LB_Partition( d_zz, 
					&d_changes,  
					&d_num_gid_entries,
					&d_num_lid_entries,
					&d_num_import,    
					&d_import_global_ids,
					&d_import_local_ids, 
					&d_import_procs,    
					&d_import_to_part,   
					&d_num_export,      
					&d_export_global_ids,
					&d_export_local_ids, 
					&d_export_procs,    
					&d_export_to_part );
#else
    Zoltan_LB_Partition( d_zz, 
			 &d_changes,  
			 &d_num_gid_entries,
			 &d_num_lid_entries,
			 &d_num_import,    
			 &d_import_global_ids,
			 &d_import_local_ids, 
			 &d_import_procs,    
			 &d_import_to_part,   
			 &d_num_export,      
			 &d_export_global_ids,
			 &d_export_local_ids, 
			 &d_export_procs,    
			 &d_export_to_part );
#endif
    testInvariant( zoltan_error == ZOLTAN_OK );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the destination process for a point given its coordinates.
 *
 * \param coords Point coordinates. The dimension of the point must be equal
 * to or less than one.
 *
 * \return The RCB destination proc for the point.
 */
template<class Mesh>
int RCB<Mesh>::getPointDestinationProc( Teuchos::Array<double> coords ) const
{
    testPrecondition( 0 <= coords.size() && coords.size() <= 3 );
    testPrecondition( d_dimension == (int) coords.size() );

    int proc = 0;
    rememberValue( int zoltan_error );
#if HAVE_DTK_DBC
    zoltan_error = Zoltan_LB_Point_Assign( d_zz, &coords[0], &proc );
#else
    Zoltan_LB_Point_Assign( d_zz, &coords[0], &proc );
#endif
    testInvariant( zoltan_error == ZOLTAN_OK );

    return proc;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the destination processes for a bounding box. This includes all
 * process domains that the box intersects.
 *
 * \param box The bounding box to get the destinations for.
 *
 * \return The RCB destination procs for the box
 */
template<class Mesh>
Teuchos::Array<int>
RCB<Mesh>::getBoxDestinationProcs( const BoundingBox& box ) const
{
    Teuchos::Tuple<double,6> box_bounds = box.getBounds();

    int num_procs = 0;
    Teuchos::Array<int> procs( d_comm->getSize() );

    rememberValue( int zoltan_error );
#if HAVE_DTK_DBC
    zoltan_error = Zoltan_LB_Box_Assign( d_zz, 
					 box[0], box[1], box[2],
					 box[3], box[4], box[5], 
					 &procs[0], &num_procs );
#else
    Zoltan_LB_Box_Assign( d_zz, 
			  box[0], box[1], box[2],
			  box[3], box[4], box[5], 
			  &procs[0], &num_procs );
#endif
    testInvariant( zoltan_error == ZOLTAN_OK );

    procs.resize( num_procs );

    return procs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the number of vertices.
 */
template<class Mesh>
int RCB<Mesh>::getNumberOfObjects( void *data, int *ierr )
{
    RCP_MeshManager mesh_manager = *static_cast<RCP_MeshManager*>( data );
    int num_vertices = 0;

    // We'll only count vertices if the mesh manager is not null.
    if ( !mesh_manager.is_null() )
    {
	int num_blocks = mesh_manager->getNumBlocks();
	Teuchos::ArrayView<short int>::const_iterator active_iterator;
	for ( int i = 0; i < num_blocks; ++i )
	{
	    for ( active_iterator = mesh_manager->getActiveVertices(i).begin();
		  active_iterator != mesh_manager->getActiveVertices(i).end();
		  ++active_iterator )
	    {
		if ( *active_iterator )
		{
		    ++num_vertices;
		}
	    }
	}
    }

    *ierr = ZOLTAN_OK;
    return num_vertices;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the local and global vertex ID's.
 */
template<class Mesh>
void RCB<Mesh>::getObjectList( 
    void *data, int sizeGID, int sizeLID,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int wgt_dim, float *obj_wgts, int *ierr )
{
    RCP_MeshManager mesh_manager = *static_cast<RCP_MeshManager*>( data );

    // We'll only build the geometry list is the mesh manager is not null.
    if ( !mesh_manager.is_null() )
    {
	// Note here that the local ID is being set the the vertex array index.
	Teuchos::ArrayRCP<short int>::const_iterator active_iterator;
	typename MT::const_vertex_iterator gid_iterator;
	zoltan_id_type i = 0;
	zoltan_id_type j = 0;
	BlockIterator block_iterator;
	for ( block_iterator = mesh_manager->blocksBegin();
	      block_iterator != mesh_manager->blocksEnd();
	      ++block_iterator )
	{
	    int block_id = std::distance( mesh_manager->blocksBegin(),
					  block_iterator );

	    for ( gid_iterator = MT::verticesBegin( *(*block_iterator) ),
	       active_iterator = 
				 mesh_manager->getActiveVertices( block_id ).begin();
		  gid_iterator != MT::verticesEnd( *(*block_iterator) );
		  ++gid_iterator, ++active_iterator )
	    {
		if ( *active_iterator )
		{
		    globalID[i] = static_cast<zoltan_id_type>( *gid_iterator );
		    localID[i] = j;
		    ++i;
		}
		++j;
	    }
	}
    }

    *ierr = ZOLTAN_OK;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the dimension of the vertices.
 */
template<class Mesh>
int RCB<Mesh>::getNumGeometry( void *data, int *ierr )
{
    int dimension = *static_cast<int*>( data );
    *ierr = ZOLTAN_OK;
    return dimension;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the vertex coordinates.
 */
template<class Mesh>
void RCB<Mesh>::getGeometryList(
    void *data, int sizeGID, int sizeLID,
    int num_obj,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int num_dim, double *geom_vec, int *ierr )
{
    RCP_MeshManager mesh_manager = *static_cast<RCP_MeshManager*>( data );

    // We will only supply vertex coordinates when the mesh exists.
    if ( !mesh_manager.is_null() )
    {
	// Get the number of active vertices.
	int num_active_vertices = 0;
	int num_blocks = mesh_manager->getNumBlocks();
	Teuchos::ArrayView<short int>::const_iterator active_iterator;
	for ( int i = 0; i < num_blocks; ++i )
	{
	    for ( active_iterator = mesh_manager->getActiveVertices(i).begin();
		  active_iterator != mesh_manager->getActiveVertices(i).end();
		  ++active_iterator )
	    {
		if ( *active_iterator )
		{
		    ++num_active_vertices;
		}
	    }
	}

	// Check Zoltan for consistency.
	int vertex_dim = mesh_manager->dim();
	testInvariant( sizeGID == 1 );
	testInvariant( sizeLID == 1 );
	testInvariant( num_dim == (int) vertex_dim );
	testInvariant( num_obj == (int) num_active_vertices );

	if ( sizeGID != 1 || sizeLID != 1 || 
	     num_dim != (int) vertex_dim || num_obj != (int) num_active_vertices )
	{
	    *ierr = ZOLTAN_FATAL;
	    return;
	}
    
	// Zoltan needs interleaved coordinates.
	int n = 0;
	Teuchos::ArrayRCP<const double> mesh_coords;
	GlobalOrdinal num_vertices;
	BlockIterator block_iterator;
	for ( block_iterator = mesh_manager->blocksBegin();
	      block_iterator != mesh_manager->blocksEnd();
	      ++block_iterator )
	{
	    int block_id = std::distance( mesh_manager->blocksBegin(),
					  block_iterator );
	    Teuchos::ArrayView<short int> active_vertices =
		mesh_manager->getActiveVertices( block_id );

	    mesh_coords = MeshTools<Mesh>::coordsView( *(*block_iterator) );
	    num_vertices = std::distance( MT::verticesBegin( *(*block_iterator) ),
					  MT::verticesEnd( *(*block_iterator) ) );
	    for ( GlobalOrdinal i = 0; i < num_vertices; ++i )
	    {
		if ( active_vertices[i] )
		{
		    for ( int d = 0; d < vertex_dim; ++d )
		    {
			geom_vec[ vertex_dim*n + d ] = 
			    mesh_coords[ d*num_vertices + i ];
		    }
		    ++n;
		}
	    }
	}
    }
	  
    *ierr = ZOLTAN_OK;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_RCB_def.hpp
//---------------------------------------------------------------------------//

