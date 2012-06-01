//---------------------------------------------------------------------------//
/*!
 * \file DTK_RCB.hpp
 * \author Stuart R. Slattery
 * \brief Wrapper definition for Zoltan recursive coordinate bisectioning.
 */
//---------------------------------------------------------------------------//

#include "DTK_RCB.hpp"
#include <DTK_Exception.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
RCB::RCB( const MPI_Comm& comm )
    : d_zz( Zoltan_Create( comm ) )
{
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
    Zoltan_Set_Num_Obj_Fn( d_zz, getNumberOfObjects, &myMesh);
    Zoltan_Set_Obj_List_Fn( d_zz, getObjectList, &myMesh);
    Zoltan_Set_Num_Geom_Fn( d_zz, getNumGeometry, &myMesh);
    Zoltan_Set_Geom_Multi_Fn( d_zz, getGeometryList, &myMesh);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
RCB::~RCB()
{
    // Tell zoltan to cleanup.
    Zoltan_LB_Free_Part ( &d_importGlobalGids, &d_importLocalGids, 
			  &d_importProcs, &d_importToPart );
    Zoltan_LB_Free_Part( &d_exportGlobalGids, &d_exportLocalGids, 
			 &d_exportProcs, &d_exportToPart );
    Zoltan_Destroy( &d_zz );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute RCB partitioning.
 */
void RCB::partition()
{
    int zoltan_error;
    zoltan_error = Zoltan_LB_Partition( d_zz, 
					&d_changes,  
					&d_numGidEntries,
					&d_numLidEntries,
					&d_numImport,    
					&d_importGlobalGids,
					&d_importLocalGids, 
					&d_importProcs,    
					&d_importToPart,   
					&d_numExport,      
					&d_exportGlobalGids,
					&d_exportLocalGids, 
					&d_exportProcs,    
					&d_exportToPart ); 

    testInvariant( ZOLTAN_OK == zoltan_error, 
		   "Zoltan error creating RCB partitioning" );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the number of nodes.
 */
static int RCB::getNumberOfObjects( void *data, int *ierr )
{

}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the local and global node ID's.
 */
static void RCB::getObjectList( void *data, int sizeGID, int sizeLID,
				ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
				int wgt_dim, float *obj_wgts, int *ierr )
{

}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the dimension of the nodes.
 */
static int RCB::getNumGeometry( void *data, int *ierr )
{

}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the node coordinates.
 */
static void RCB::getGeometryList( void *data, int sizeGID, int sizeLID,
				  int num_obj,
				  ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
				  int num_dim, double *geom_vec, int *ierr )
{

}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_RCB.hpp
//---------------------------------------------------------------------------//

