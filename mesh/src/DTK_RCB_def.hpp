//---------------------------------------------------------------------------//
/*!
 * \file DTK_RCB_def.hpp
 * \author Stuart R. Slattery
 * \brief Wrapper definition for Zoltan recursive coordinate bisectioning.
 */
//---------------------------------------------------------------------------//

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
    : d_mesh_data( mesh, active_nodes )
    , d_part_boxes( comm->getSize() )
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

    // Get the bounding boxes on all processes.
    Teuchos::Array<BoundingBox>::iterator box_iterator;
    int i = 0;
    for ( box_iterator = d_part_boxes.begin();
	  box_iterator != d_part_boxes.end();
	  ++box_iterator, ++i )
    {
	*box_iterator = getPartBoundingBox( i );
    }

    // We should build an octree or other logarthmic search structure here out
    // of the bounding boxes. That way getDestinationProc() will avoid linear
    // searches (and ultimately quadratic performance).
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the destination process for a point.
 */
template<class Mesh>
int RCB<Mesh>::getDestinationProc( double coords[3] ) const
{
    // Do a linear search through the bounding boxes for now. This really
    // needs to be logarithmic as we are checking every mesh node with this
    // during rendezvous construction and every target node during mapping.
    Teuchos::Array<BoundingBox>::const_iterator box_iterator;
    int i = 0;
    for ( box_iterator = d_part_boxes.begin();
	  box_iterator != d_part_boxes.end();
	  ++box_iterator, ++i )
    {
	if ( box_iterator->pointInBox( coords ) )
	{
	    return i;
	}
    }

    // We didn't find the point in the RCB partitioning so throw a point not
    // found exception.
    throw PointNotFound( "Did not find point in the RCB decomposition." );

    return 0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the bounding box for a partition.
 */
template<class Mesh>
BoundingBox RCB<Mesh>::getPartBoundingBox( const int part ) const
{
    double x_min, y_min, z_min, x_max, y_max, z_max;
    int dim;
    int zoltan_error;
    zoltan_error = Zoltan_RCB_Box( d_zz, part, &dim,
				   &x_min, &y_min, &z_min,
				   &x_max, &y_max, &z_max );

    testInvariant( ZOLTAN_OK == zoltan_error, 
		   "Zoltan error getting partition bounding box." );

    return BoundingBox( x_min, y_min, z_min, x_max, y_max, z_max );
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
    typename MT::const_node_iterator handle_iterator;
    zoltan_id_type i = 0;
    zoltan_id_type j = 0;
    for ( handle_iterator = MT::nodesBegin( mesh_data->d_mesh ),
	  active_iterator = mesh_data->d_active_nodes.begin();
	  handle_iterator != MT::nodesEnd( mesh_data->d_mesh );
	  ++handle_iterator, ++active_iterator )
    {
	if ( *active_iterator )
	{
	    globalID[i] = static_cast<zoltan_id_type>( *handle_iterator );
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
    typename MT::const_coordinate_iterator mesh_coords =
	MT::coordsBegin( mesh_data->d_mesh );
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

