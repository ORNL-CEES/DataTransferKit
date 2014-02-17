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
 * \file DTK_GeometryRCB.hpp
 * \author Stuart R. Slattery
 * \brief Wrapper declaration for Zoltan recursive coordinate bisectioning.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_GEOMETRYRCB_HPP
#define DTK_GEOMETRYRCB_HPP

#include <set>

#include "DTK_Partitioner.hpp"
#include "DTK_BoundingBox.hpp"
#include "DTK_GeometryTraits.hpp"
#include "DTK_GeometryManager.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <zoltan.h>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \class GeometryRCB
 * \brief Recursive Coordinate Bisectioning partitioner for geometry
 */
//---------------------------------------------------------------------------//
template<class Geometry, class GlobalOrdinal>
class GeometryRCB : public Partitioner
{
  public:
    
    //@{
    //! Typedefs.
    typedef Geometry                                     geometry_type;
    typedef GlobalOrdinal                                global_ordinal_type;
    typedef GeometryTraits<Geometry>                     GT;
    typedef GeometryManager<Geometry,GlobalOrdinal>      GeometryManagerType;
    typedef Teuchos::RCP<GeometryManagerType>            RCP_GeometryManager;
    typedef Teuchos::Comm<int>                           CommType;
    typedef Teuchos::RCP<const CommType>                 RCP_Comm;
    typedef ZOLTAN_ID_TYPE                               zoltan_id_type;
    typedef ZOLTAN_ID_PTR                                zoltan_id_ptr;
    //@}

    // Constructor.
    GeometryRCB( const RCP_Comm& comm, 
		 const RCP_GeometryManager& geometry_manager, 
		 const int dimension );

    // Destructor.
    ~GeometryRCB();

    // Compute GeometryRCB partitioning of the geometry.
    void partition( const BoundingBox& local_box );

    // Given a range of local input point ids get their destination procs.
    Teuchos::Array<int> getInputPointDestinationProcs(
	const int lid_begin, const int num_points );

    // Get the destination process for a point given its coordinates.
    int getPointDestinationProc( Teuchos::ArrayView<double> coords ) const;

    // Get the destination processes for a bounding box.
    Teuchos::Array<int> getBoxDestinationProcs( const BoundingBox& box ) const;

    //! Get the number of imported vertices.
    int getNumImport() const
    { return d_num_import; }

    //! Get the global import vertex ID's.
    Teuchos::ArrayView<zoltan_id_type> getImportGlobalIds() const
    { return Teuchos::ArrayView<zoltan_id_type>( d_import_global_ids, 
						 d_num_import ); }

    //! Get the local import vertex ID's.
    Teuchos::ArrayView<zoltan_id_type> getImportLocalIds() const
    { return Teuchos::ArrayView<zoltan_id_type>( d_import_local_ids, 
						 d_num_import ); }

    //! Get the process rank source for imported vertices.
    Teuchos::ArrayView<int> getImportProcs() const
    { return Teuchos::ArrayView<int>( d_import_procs, d_num_import ); }

    //! Get the new partition for imported vertices.
    Teuchos::ArrayView<int> getImportParts() const
    { return Teuchos::ArrayView<int>( d_import_to_part, d_num_import ); }

    //! Get the number of exported vertices.
    int getNumExport() const
    { return d_num_export; }

    //! Get the global export vertex ID's.
    Teuchos::ArrayView<zoltan_id_type> getExportGlobalIds() const
    { return Teuchos::ArrayView<zoltan_id_type>( d_export_global_ids, 
						 d_num_export ); }

    //! Get the local export vertex ID's.
    Teuchos::ArrayView<zoltan_id_type> getExportLocalIds() const
    { return Teuchos::ArrayView<zoltan_id_type>( d_export_local_ids, 
						 d_num_export ); }

    //! Get the process rank target for exported vertices.
    Teuchos::ArrayView<int> getExportProcs() const
    { return Teuchos::ArrayView<int>( d_export_procs, d_num_export ); }

    //! Get the new partition for exported vertices.
    Teuchos::ArrayView<int> getExportParts() const
    { return Teuchos::ArrayView<int>( d_export_to_part, d_num_export ); }

  private:

    // Zoltan callback for getting the number of vertices.
    static int getNumberOfObjects( void *data, int *ierr );

    // Zoltan callback for getting the local and global vertex ID's.
    static void getObjectList( void *data, int sizeGID, int sizeLID,
			       ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			       int wgt_dim, float *obj_wgts, int *ierr );

    // Zoltan callback for getting the dimension of the vertices.
    static int getNumGeometry( void *data, int *ierr );

    // Zoltan callback for getting the vertex coordinates.
    static void getGeometryList( void *data, int sizeGID, int sizeLID,
				 int num_obj,
				 ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
				 int num_dim, double *geom_vec, int *ierr );

  private:

    // The communicator over which GeometryRCB is performed.
    RCP_Comm d_comm;

    // The geometry being partitioned.
    RCP_GeometryManager d_geometry_manager;

    // The dimension of the GeometryRCB space.
    int d_dimension;

    // Bounding boxes for source processes that neighbor this process.
    Teuchos::Array<BoundingBox> d_rcb_boxes;

    // Ranks for source processes that neighbor this process.
    Teuchos::Array<int> d_box_ranks;

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

#include "DTK_GeometryRCB_def.hpp"

//---------------------------------------------------------------------------//

#endif // end DTK_GEOMETRYRCB_HPP

//---------------------------------------------------------------------------//
// end DTK_GeometryRCB.hpp
//---------------------------------------------------------------------------//
