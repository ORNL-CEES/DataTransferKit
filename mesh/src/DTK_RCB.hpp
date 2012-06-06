//---------------------------------------------------------------------------//
/*!
 * \file DTK_RCB.hpp
 * \author Stuart R. Slattery
 * \brief Wrapper declaration for Zoltan recursive coordinate bisectioning.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_RCB_HPP
#define DTK_RCB_HPP

#include <vector>

#include "DTK_BoundingBox.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>

#include <zoltan.h>

namespace DataTransferKit
{

template<typename Mesh>
class RCB
{

  public:
    
    //@{
    //! Typedefs.
    typedef Mesh                                        mesh_type;
    typedef Teuchos::Comm<int>                          CommType;
    typedef Teuchos::RCP<const CommType>                RCP_Comm;
    //@}

    // Constructor.
    RCB( const Mesh& mesh, const RCP_Comm& comm );

    // Destructor.
    ~RCB();

    // Compute RCB partitioning of the node field.
    void partition();

    // Get the destination process for a point.
    int getDestinationProc( double coords[3] ) const;

    //! Get the number of imported nodes.
    int getNumImport() const
    { return d_num_import; }

    //! Get the global import node ID's.
    Teuchos::ArrayView<unsigned int> getImportGlobalIds() const
    { return Teuchos::ArrayView<unsigned int>( d_import_global_ids, 
					       d_num_import ); }

    //! Get the local import node ID's.
    Teuchos::ArrayView<unsigned int> getImportLocalIds() const
    { return Teuchos::ArrayView<unsigned int>( d_import_local_ids, 
					       d_num_import ); }

    //! Get the process rank source for imported nodes.
    Teuchos::ArrayView<int> getImportProcs() const
    { return Teuchos::ArrayView<int>( d_import_procs, d_num_import ); }

    //! Get the new partition for imported nodes.
    Teuchos::ArrayView<int> getImportParts() const
    { return Teuchos::ArrayView<int>( d_import_to_part, d_num_import ); }

    //! Get the number of exported nodes.
    int getNumExport() const
    { return d_num_export; }

    //! Get the global export node ID's.
    Teuchos::ArrayView<unsigned int> getExportGlobalIds() const
    { return Teuchos::ArrayView<unsigned int>( d_export_global_ids, 
					       d_num_export ); }

    //! Get the local export node ID's.
    Teuchos::ArrayView<unsigned int> getExportLocalIds() const
    { return Teuchos::ArrayView<unsigned int>( d_export_local_ids, 
					       d_num_export ); }

    //! Get the process rank source for exported nodes.
    Teuchos::ArrayView<int> getExportProcs() const
    { return Teuchos::ArrayView<int>( d_export_procs, d_num_export ); }

    //! Get the new partition for exported nodes.
    Teuchos::ArrayView<int> getExportParts() const
    { return Teuchos::ArrayView<int>( d_export_to_part, d_num_export ); }

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

    // Get the bounding box for a partition.
    BoundingBox getPartBoundingBox( const int part ) const;

  private:

    // The mesh we are partitioning.
    Mesh d_mesh;

    // Partition bounding boxes.
    std::vector<BoundingBox> d_part_boxes;

    // Zoltan struct.
    Zoltan_Struct *d_zz;

    // 1 if partitioning was changed, 0 otherwise.
    int d_changes;

    // Number of integers used for a global ID.
    int d_num_gid_entries;

    // Number of integers used for a local ID.
    int d_num_lid_entries;

    // Number of vertices to be sent to me.
    int d_num_import;

    // Global IDs of vertices to be sent to me.
    ZOLTAN_ID_PTR d_import_global_ids;

    // Local IDs of vertices to be sent to me.
    ZOLTAN_ID_PTR d_import_local_ids;

    // Process rank for source of each incoming vertex.
    int *d_import_procs;

    // New partition for each incoming vertex.
    int *d_import_to_part;

    // Number of vertices I must send to other processes.
    int d_num_export;

    // Global IDs of the vertices I must send.
    ZOLTAN_ID_PTR d_export_global_ids;

    // Local IDs of the vertices I must send.
    ZOLTAN_ID_PTR d_export_local_ids;

    // Process to which I send each of the vertices.
    int *d_export_procs;
    
    // Partition to which each vertex will belong.
    int *d_export_to_part;
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
