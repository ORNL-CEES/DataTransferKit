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
#include <DTK_Exception.hpp>

#include <mpi.h>

#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_CommHelpers.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Mesh>
RCB<Mesh>::RCB( const Mesh& mesh, const Teuchos::ArrayRCP<int>& active_nodes,
		const RCP_Comm& comm )
    : d_comm( comm )
    , d_mesh_data( mesh, active_nodes )
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
    Zoltan_Set_Param( d_zz, "RCB_SET_DIRECTIONS", "1" );

    // Register query functions.
    Zoltan_Set_Num_Obj_Fn( d_zz, getNumberOfObjects, &d_mesh_data );
    Zoltan_Set_Obj_List_Fn( d_zz, getObjectList, &d_mesh_data );
    Zoltan_Set_Num_Geom_Fn( d_zz, getNumGeometry, &d_mesh_data );
    Zoltan_Set_Geom_Multi_Fn( d_zz, getGeometryList, &d_mesh_data );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Mesh>
RCB<Mesh>::~RCB()
{
    // Zoltan cleanup.
    Zoltan_LB_Free_Part ( &d_import_global_ids, &d_import_local_ids, 
			  &d_import_procs, &d_import_to_part );
    Zoltan_LB_Free_Part( &d_export_global_ids, &d_export_local_ids, 
			 &d_export_procs, &d_export_to_part );
    Zoltan_Destroy( &d_zz );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute RCB partitioning of the node field.
 */
template<class Mesh>
void RCB<Mesh>::partition()
{
    // Run zoltan partitioning.
    int zoltan_error;
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

    testInvariant( ZOLTAN_OK == zoltan_error, 
		   "Zoltan error creating RCB partitioning" );

    // Get the global partitioning information.
    getPartitioning();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the destination process for a point.
 */
template<class Mesh>
int RCB<Mesh>::getDestinationProc( double coords[3] ) const
{
    int x_idx = std::distance( d_x_edges.begin(),
			       std::upper_bound( d_x_edges.begin(),
						 d_x_edges.end(),
						 coords[0] ) );

    int y_idx = std::distance( d_y_edges.begin(),
			       std::upper_bound( d_y_edges.begin(),
						 d_y_edges.end(),
						 coords[1] ) );

    int z_idx = std::distance( d_z_edges.begin(),
			       std::upper_bound( d_z_edges.begin(),
						 d_z_edges.end(),
						 coords[2] ) );

    if ( x_idx == 0 || y_idx == 0 || z_idx == 0 ||
	 x_idx > (int) d_x_edges.size() || y_idx > (int) d_y_edges.size() ||
	 z_idx > (int) d_z_edges.size() )
    {
	throw PointNotFound( "Point outside RCB decomposition." );
    }

    return (x_idx-1) + (y_idx-1)*(d_x_edges.size()-1) + 
	(z_idx-1)*(d_y_edges.size()-1)*(d_x_edges.size()-1);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the global partitioning information
 */
template<class Mesh>
void RCB<Mesh>::getPartitioning()
{
    double x_min, y_min, z_min, x_max, y_max, z_max;
    int dim;
    int zoltan_error;

    for ( int i = 0; i < d_comm->getSize(); ++i )
    {
	zoltan_error = Zoltan_RCB_Box( d_zz, i, &dim,
				       &x_min, &y_min, &z_min,
				       &x_max, &y_max, &z_max );

	testInvariant( ZOLTAN_OK == zoltan_error, 
		       "Zoltan error getting partition bounding box." );

	d_x_edges.insert( x_min );
	d_x_edges.insert( x_max );

	d_y_edges.insert( y_min );
	d_y_edges.insert( y_max );

	d_z_edges.insert( z_min );
	d_z_edges.insert( z_max );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the number of nodes.
 */
template<class Mesh>
int RCB<Mesh>::getNumberOfObjects( void *data, int *ierr )
{
    MeshData *mesh_data = static_cast<MeshData*>( data );
    int num_nodes = 0;
    Teuchos::ArrayRCP<int>::const_iterator active_iterator;
    for ( active_iterator = mesh_data->d_active_nodes.begin();
	  active_iterator != mesh_data->d_active_nodes.end();
	  ++active_iterator )
    {
	if ( *active_iterator )
	{
	    ++num_nodes;
	}
    }

    *ierr = ZOLTAN_OK;
    return num_nodes;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the local and global node ID's.
 */
template<class Mesh>
void RCB<Mesh>::getObjectList( 
    void *data, int sizeGID, int sizeLID,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int wgt_dim, float *obj_wgts, int *ierr )
{
    MeshData *mesh_data = static_cast<MeshData*>( data );
    *ierr = ZOLTAN_OK;

    // Note here that the local ID is being set the the node array index.
    Teuchos::ArrayRCP<int>::const_iterator active_iterator;
    typename MT::const_node_iterator gid_iterator;
    zoltan_id_type i = 0;
    zoltan_id_type j = 0;
    for ( gid_iterator = MT::nodesBegin( mesh_data->d_mesh ),
       active_iterator = mesh_data->d_active_nodes.begin();
	  gid_iterator != MT::nodesEnd( mesh_data->d_mesh );
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

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the dimension of the nodes.
 */
template<class Mesh>
int RCB<Mesh>::getNumGeometry( void *data, int *ierr )
{
    MeshData *mesh_data = static_cast<MeshData*>( data );
    *ierr = ZOLTAN_OK;
    return MT::nodeDim( mesh_data->d_mesh );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the node coordinates.
 */
template<class Mesh>
void RCB<Mesh>::getGeometryList(
    void *data, int sizeGID, int sizeLID,
    int num_obj,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int num_dim, double *geom_vec, int *ierr )
{
    MeshData *mesh_data = static_cast<MeshData*>( data );

    // Get the number of active nodes.
    std::size_t num_nodes = 0;
    Teuchos::ArrayRCP<int>::const_iterator active_iterator;
    for ( active_iterator = mesh_data->d_active_nodes.begin();
	  active_iterator != mesh_data->d_active_nodes.end();
	  ++active_iterator )
    {
	if ( *active_iterator )
	{
	    ++num_nodes;
	}
    }

    // Check Zoltan for consistency.
    std::size_t node_dim = MT::nodeDim( mesh_data->d_mesh );
    testInvariant( sizeGID == 1, "Zoltan global ID size != 1." );
    testInvariant( sizeLID == 1, "Zoltan local ID size != 1." );
    testInvariant( num_dim == (int) node_dim, "Zoltan dimension != 3." );
    testInvariant( num_obj == (int) num_nodes, 
		   "Zoltan number of nodes != mesh number of nodes." );

    if ( sizeGID != 1 || sizeLID != 1 || 
	 num_dim != (int) node_dim || num_obj != (int) num_nodes )
    {
	*ierr = ZOLTAN_FATAL;
	return;
    }
    
    // Zoltan needs interleaved coordinates.
    Teuchos::ArrayRCP<const double> mesh_coords = 
	MeshTools<Mesh>::coordsView( mesh_data->d_mesh );
    for ( std::size_t n = 0; n < num_nodes; ++n )
    {
	for ( std::size_t d = 0; d < node_dim; ++d )
	{
	    geom_vec[ node_dim*n + d ] = mesh_coords[ d*num_nodes + n ];
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_RCB_def.hpp
//---------------------------------------------------------------------------//

