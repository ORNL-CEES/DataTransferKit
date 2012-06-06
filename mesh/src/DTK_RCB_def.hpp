//---------------------------------------------------------------------------//
/*!
 * \file DTK_RCB_def.hpp
 * \author Stuart R. Slattery
 * \brief Wrapper definition for Zoltan recursive coordinate bisectioning.
 */
//---------------------------------------------------------------------------//

#include <DTK_MeshTraits.hpp>
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
template<typename MeshType>
RCB<MeshType>::RCB( const MeshType& mesh, const RCP_Comm& comm )
    : d_mesh( mesh )
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
    Zoltan_Set_Num_Obj_Fn( d_zz, getNumberOfObjects, &d_mesh );
    Zoltan_Set_Obj_List_Fn( d_zz, getObjectList, &d_mesh );
    Zoltan_Set_Num_Geom_Fn( d_zz, getNumGeometry, &d_mesh );
    Zoltan_Set_Geom_Multi_Fn( d_zz, getGeometryList, &d_mesh );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<typename MeshType>
RCB<MeshType>::~RCB()
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
template<typename MeshType>
void RCB<MeshType>::partition()
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
    std::vector<BoundingBox>::iterator box_iterator;
    int i = 0;
    for ( box_iterator = d_part_boxes.begin();
	  box_iterator != d_part_boxes.end();
	  ++box_iterator, ++i )
    {
	*box_iterator = getPartBoundingBox( i );
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the bounding box for a partition.
 */
template<typename MeshType>
BoundingBox RCB<MeshType>::getPartBoundingBox( const int part ) const
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
 * \brief Get the destination process for a point.
 */
template<typename MeshType>
int RCB<MeshType>::getDestinationProc( double coords[3] ) const
{
    // Do a linear search through the bounding boxes for now.
    std::vector<BoundingBox>::const_iterator box_iterator;
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
 * \brief Zoltan callback for getting the number of nodes.
 */
template<typename MeshType>
int RCB<MeshType>::getNumberOfObjects( void *data, int *ierr )
{
    MeshType *mesh = (MeshType*) data;
    *ierr = ZOLTAN_OK;
    return std::distance( MeshTraits<MeshType>::nodesBegin( *mesh ),
			  MeshTraits<MeshType>::nodesEnd( *mesh ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the local and global node ID's.
 */
template<typename MeshType>
void RCB<MeshType>::getObjectList( 
    void *data, int sizeGID, int sizeLID,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int wgt_dim, float *obj_wgts, int *ierr )
{
    MeshType *mesh = (MeshType*) data;
    *ierr = ZOLTAN_OK;

    typename MeshTraits<MeshType>::const_handle_iterator handle_iterator;
    int i = 0;
    for ( handle_iterator = MeshTraits<MeshType>::nodesBegin( *mesh );
	  handle_iterator != MeshTraits<MeshType>::nodesEnd( *mesh );
	  ++handle_iterator, ++i )
    {
	globalID[i] = (ZOLTAN_ID_TYPE) *handle_iterator;
	localID[i] = i;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the dimension of the nodes.
 */
template<typename MeshType>
int RCB<MeshType>::getNumGeometry( void *data, int *ierr )
{
    *ierr = ZOLTAN_OK;
    return 3;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the node coordinates.
 */
template<typename MeshType>
void RCB<MeshType>::getGeometryList(
    void *data, int sizeGID, int sizeLID,
    int num_obj,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int num_dim, double *geom_vec, int *ierr )
{
    MeshType *mesh = (MeshType*) data;
    int num_nodes = std::distance( MeshTraits<MeshType>::nodesBegin( *mesh ),
				   MeshTraits<MeshType>::nodesEnd( *mesh ) );

    testInvariant( sizeGID == 1, "Zoltan global ID size != 1." );
    testInvariant( sizeLID == 1, "Zoltan local ID size != 1." );
    testInvariant( num_dim == 3, "Zoltan dimension != 3." );
    testInvariant( num_obj == num_nodes, 
		   "Zoltan number of nodes != mesh number of nodes." );

    if ( sizeGID != 1 || sizeLID != 1 || num_dim != 3 || num_obj != num_nodes )
    {
	*ierr = ZOLTAN_FATAL;
	return;
    }
    
    typename MeshTraits<MeshType>::const_coordinate_iterator coord_iterator;
    int i = 0;
    for ( coord_iterator = MeshTraits<MeshType>::coordsBegin( *mesh );
	  coord_iterator != MeshTraits<MeshType>::coordsEnd( *mesh );
	  ++coord_iterator )
    {
	geom_vec[i] = *coord_iterator;
	++i;
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_RCB_def.hpp
//---------------------------------------------------------------------------//

