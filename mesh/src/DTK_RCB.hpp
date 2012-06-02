//---------------------------------------------------------------------------//
/*!
 * \file DTK_RCB.hpp
 * \author Stuart R. Slattery
 * \brief Wrapper declaration for Zoltan recursive coordinate bisectioning.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_RCB_HPP
#define DTK_RCB_HPP

#include <mpi.h>

#include <zoltan.h>

namespace DataTransferKit
{

template<typename NodeField>
class RCB
{

  public:

    // Constructor.
    RCB( const NodeField& node_field, const MPI_Comm& comm );

    // Destructor.
    ~RCB();

    // Compute RCB partitioning of the node field.
    void partition();

  private:

    // Zoltan callback for getting the number of nodes.
    static int getNumberOfObjects( void *data, int *ierr );

    // Zoltan callback for getting the local and global node ID's.
    static void getObjectList( void *data, int sizeGID, int sizeLID,
			       ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			       int wgt_dim, float *obj_wgts, int *ierr );

    // Zoltan callback for getting the dimension of the nodes.
    static int getNumGeometry( void *data, int *ierr );

    // Zoltan callback for getting the node coordinates.
    static void getGeometryList( void *data, int sizeGID, int sizeLID,
				 int num_obj,
				 ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
				 int num_dim, double *geom_vec, int *ierr );

  private:

    // The node field we are partitioning.
    NodeField d_node_field;

    // Zoltan struct.
    Zoltan_Struct *d_zz;

    // 1 if partitioning was changed, 0 otherwise.
    int d_changes;

    // Number of integers used for a global ID.
    int d_numGidEntries;

    // Number of integers used for a local ID.
    int d_numLidEntries;

    // Number of vertices to be sent to me.
    int d_numImport;

    // Global IDs of vertices to be sent to me.
    ZOLTAN_ID_PTR d_importGlobalGids;

    // Local IDs of vertices to be sent to me.
    ZOLTAN_ID_PTR d_importLocalGids;

    // Process rank for source of each incoming vertex.
    int *d_importProcs;

    // New partition for each incoming vertex.
    int *d_importToPart;

    // Number of vertices I must send to other processes.
    int d_numExport;

    // Global IDs of the vertices I must send.
    ZOLTAN_ID_PTR d_exportGlobalGids;

    // Local IDs of the vertices I must send.
    ZOLTAN_ID_PTR d_exportLocalGids;

    // Process to which I send each of the vertices.
    int *d_exportProcs;
    
    // Partition to which each vertex will belong.
    int *d_exportToPart;
};

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "DTK_RCB_def.hpp"

#endif // end DTK_RCB_HPP

//---------------------------------------------------------------------------//
// end DTK_RCB.hpp
//---------------------------------------------------------------------------//
