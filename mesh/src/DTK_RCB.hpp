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
 * \file DTK_RCB.hpp
 * \author Stuart R. Slattery
 * \brief Wrapper declaration for Zoltan recursive coordinate bisectioning.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_RCB_HPP
#define DTK_RCB_HPP

#include <set>

#include "DTK_BoundingBox.hpp"
#include "DTK_MeshTraits.hpp"
#include "DTK_MeshManager.hpp"

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
 * \class RCB
 * \brief Recursive Coordinate Bisectioning partitioner.
 */
//---------------------------------------------------------------------------//
template<class Mesh>
class RCB
{

  public:
    
    //@{
    //! Typedefs.
    typedef Mesh                                        mesh_type;
    typedef MeshTraits<Mesh>                            MT;
    typedef typename MT::global_ordinal_type            GlobalOrdinal;
    typedef Teuchos::RCP< MeshManager<Mesh> >           RCP_MeshManager;
    typedef typename MeshManager<Mesh>::BlockIterator   BlockIterator;          
    typedef Teuchos::Comm<int>                          CommType;
    typedef Teuchos::RCP<const CommType>                RCP_Comm;
    typedef ZOLTAN_ID_TYPE                              zoltan_id_type;
    typedef ZOLTAN_ID_PTR                               zoltan_id_ptr;
    //@}

    // Constructor.
    RCB( const RCP_MeshManager& mesh_manager );

    // Destructor.
    ~RCB();

    // Compute RCB partitioning of the mesh.
    void partition();

    // Get the destination process for a node given its coordinates.
    int getDestinationProc( const Teuchos::Array<double>& coords ) const;

    //! Get the number of imported nodes.
    int getNumImport() const
    { return d_num_import; }

    //! Get the global import node ID's.
    Teuchos::ArrayView<zoltan_id_type> getImportGlobalIds() const
    { return Teuchos::ArrayView<zoltan_id_type>( d_import_global_ids, 
						 d_num_import ); }

    //! Get the local import node ID's.
    Teuchos::ArrayView<zoltan_id_type> getImportLocalIds() const
    { return Teuchos::ArrayView<zoltan_id_type>( d_import_local_ids, 
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
    Teuchos::ArrayView<zoltan_id_type> getExportGlobalIds() const
    { return Teuchos::ArrayView<zoltan_id_type>( d_export_global_ids, 
						 d_num_export ); }

    //! Get the local export node ID's.
    Teuchos::ArrayView<zoltan_id_type> getExportLocalIds() const
    { return Teuchos::ArrayView<zoltan_id_type>( d_export_local_ids, 
						 d_num_export ); }

    //! Get the process rank target for exported nodes.
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

    // Get the global partitioning information.
    void getPartitioning();

  private:

    // The mesh to be partitioned.
    RCP_MeshManager d_mesh_manager;

    // Global x edges.
    std::set<double> d_x_edges;

    // Global y edges.
    std::set<double> d_y_edges;

    // Global z edges.
    std::set<double> d_z_edges;

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
