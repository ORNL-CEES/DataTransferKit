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
    NodeField *node_field = (NodeField*) data;
    *ierr = ZOLTAN_OK;
    return FieldTraits<NodeField>::size( *node_field );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the local and global node ID's.
 */
template<typename NodeField>
void RCB<NodeField>::getObjectList( 
    void *data, int sizeGID, int sizeLID,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int wgt_dim, float *obj_wgts, int *ierr )
{
    NodeField *node_field = (NodeField*) data;
    *ierr = ZOLTAN_OK;

    typedef typename FieldTraits<NodeField>::value_type node_type;
    
    typename FieldTraits<NodeField>::const_iterator node_iterator;
    int i = 0;
    for ( node_iterator = FieldTraits<NodeField>::begin( *node_field );
	  node_iterator != FieldTraits<NodeField>::begin( *node_field );
	  ++node_iterator, ++i )
    {
	globalID[i] = 
	    (ZOLTAN_ID_TYPE) NodeTraits<node_type>::handle( *node_iterator );
	localID[i] = i;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the dimension of the nodes.
 */
template<typename NodeField>
int RCB<NodeField>::getNumGeometry( void *data, int *ierr )
{
    *ierr = ZOLTAN_OK;
    typedef typename FieldTraits<NodeField>::value_type node_type;
    return NodeTraits<node_type>::dim();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the node coordinates.
 */
template<typename NodeField>
void RCB<NodeField>::getGeometryList(
    void *data, int sizeGID, int sizeLID,
    int num_obj,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int num_dim, double *geom_vec, int *ierr )
{
    NodeField *node_field = (NodeField*) data;
    int num_nodes = FieldTraits<NodeField>::size( *node_field );

    typedef typename FieldTraits<NodeField>::value_type node_type;
    int dim = NodeTraits<node_type>::dim();

    testInvariant( sizeGID == 1, "Zoltan global ID size != 1." );
    testInvariant( sizeLID == 1, "Zoltan local ID size != 1." );
    testInvariant( num_dim == dim, "Zoltan dimension != node dimension." );
    testInvariant( num_obj == num_nodes, 
		   "Zoltan number of nodes != field size" );

    if ( sizeGID != 1 || sizeLID != 1 || num_dim != dim || num_obj != num_nodes )
    {
	*ierr = ZOLTAN_FATAL;
	return;
    }
    
    typename FieldTraits<NodeField>::const_iterator node_iterator;
    typename NodeTraits<node_type>::const_coordinate_iterator coord_iterator;
    int i = 0;
    for ( node_iterator = FieldTraits<NodeField>::begin( *node_field );
	  node_iterator != FieldTraits<NodeField>::begin( *node_field );
	  ++node_iterator)
    {
	for ( coord_iterator = 
		  NodeTraits<node_type>::coordsBegin( *node_iterator );
	      coord_iterator != 
		  NodeTraits<node_type>::coordsEnd( *node_iterator );
	      ++coord_iterator )
	{
	    geom_vec[i] = (double) *coord_iterator;
	    ++i;
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// end DTK_RCB_def.hpp
//---------------------------------------------------------------------------//

