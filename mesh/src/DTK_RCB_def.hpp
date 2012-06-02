//---------------------------------------------------------------------------//
/*!
 * \file DTK_RCB_def.hpp
 * \author Stuart R. Slattery
 * \brief Wrapper definition for Zoltan recursive coordinate bisectioning.
 */
//---------------------------------------------------------------------------//

#include <DTK_Exception.hpp>
#include <DTK_NodeTraits.hpp>
#include <DTK_FieldTraits.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename NodeField>
RCB<NodeField>::RCB( const NodeField& node_field, const MPI_Comm& comm )
    : d_node_field( node_field )
    , d_zz( Zoltan_Create( comm ) )
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
    Zoltan_Set_Num_Obj_Fn( d_zz, getNumberOfObjects, &d_node_field );
    Zoltan_Set_Obj_List_Fn( d_zz, getObjectList, &d_node_field );
    Zoltan_Set_Num_Geom_Fn( d_zz, getNumGeometry, &d_node_field );
    Zoltan_Set_Geom_Multi_Fn( d_zz, getGeometryList, &d_node_field );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<typename NodeField>
RCB<NodeField>::~RCB()
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
 * \brief Compute RCB partitioning of the node field.
 */
template<typename NodeField>
void RCB<NodeField>::partition()
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
template<typename NodeField>
int RCB<NodeField>::getNumberOfObjects( void *data, int *ierr )
{

}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the local and global node ID's.
 */
template<typename NodeField>
void RCB<NodeField>::getObjectList( void *data, int sizeGID, int sizeLID,
				    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
				    int wgt_dim, float *obj_wgts, int *ierr )
{

}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the dimension of the nodes.
 */
template<typename NodeField>
int RCB<NodeField>::getNumGeometry( void *data, int *ierr )
{

}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the node coordinates.
 */
template<typename NodeField>
void RCB<NodeField>::getGeometryList( void *data, int sizeGID, int sizeLID,
				      int num_obj,
				      ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
				      int num_dim, double *geom_vec, int *ierr )
{

}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_RCB_def.hpp
//---------------------------------------------------------------------------//

