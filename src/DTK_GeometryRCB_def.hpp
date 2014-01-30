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
 * \file DTK_GeometryRCB_def.hpp
 * \author Stuart R. Slattery
 * \brief Wrapper definition for Zoltan recursive coordinate bisectioning.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_GEOMETRYRCB_DEF_HPP
#define DTK_GEOMETRYRCB_DEF_HPP

#include <numeric>
#include <algorithm>

#include "DTK_DBC.hpp"
#include "DTK_CommIndexer.hpp"
#include "DataTransferKit_config.hpp"

#include <mpi.h>

#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_as.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param comm The communicator over which to build the GeometryRCB partitioning.
 *
 * \param geometry_manager The geometry to be partitioned with GeometryRCB. A
 * null RCP is valid here as the geometry may or may not exist on all of the
 * processes we want to repartition it to.
 *
 * \param dimension The dimension of the GeometryRCB space.
 */
template<class Geometry, class GlobalOrdinal>
GeometryRCB<Geometry,GlobalOrdinal>::GeometryRCB(
    const RCP_Comm& comm, const RCP_GeometryManager& geometry_manager, 
    const int dimension )
    : d_comm( comm )
    , d_geometry_manager( geometry_manager )
    , d_dimension( dimension )
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

    // GeometryRCB parameters.
    Zoltan_Set_Param( d_zz, "RCB_OUTPUT_LEVEL", "0" );
    Zoltan_Set_Param( d_zz, "RCB_RECTILINEAR_BLOCKS", "0" );
    Zoltan_Set_Param( d_zz, "KEEP_CUTS", "1" );
    Zoltan_Set_Param( d_zz, "AVERAGE_CUTS", "1" );
    Zoltan_Set_Param( d_zz, "RCB_LOCK_DIRECTIONS", "1" );
    Zoltan_Set_Param( d_zz, "RCB_SET_DIRECTIONS", "1" );

    // Register static functions.
    Zoltan_Set_Num_Obj_Fn( d_zz, getNumberOfObjects, &d_geometry_manager );
    Zoltan_Set_Obj_List_Fn( d_zz, getObjectList, &d_geometry_manager );
    Zoltan_Set_Num_Geom_Fn( d_zz, getNumGeometry, &d_dimension );
    Zoltan_Set_Geom_Multi_Fn( d_zz, getGeometryList, &d_geometry_manager );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor. Zoltan memory deallocation happens here and only here.
 */
template<class Geometry, class GlobalOrdinal>
GeometryRCB<Geometry,GlobalOrdinal>::~GeometryRCB()
{
    Zoltan_LB_Free_Part( &d_import_global_ids, &d_import_local_ids, 
			 &d_import_procs, &d_import_to_part );
    Zoltan_LB_Free_Part( &d_export_global_ids, &d_export_local_ids, 
			 &d_export_procs, &d_export_to_part );
    Zoltan_Destroy( &d_zz );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute GeometryRCB partitioning of the geometry.
 */
template<class Geometry, class GlobalOrdinal>
void GeometryRCB<Geometry,GlobalOrdinal>::partition()
{
    DTK_REMEMBER( int zoltan_error );
#if HAVE_DTK_DBC
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
#else
    Zoltan_LB_Partition( d_zz, 
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
#endif
    DTK_CHECK( zoltan_error == ZOLTAN_OK );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the destination process for a point given its coordinates.
 *
 * \param coords Point coordinates. The dimension of the point must be equal
 * to or less than one.
 *
 * \return The GeometryRCB destination proc for the point.
 */
template<class Geometry, class GlobalOrdinal>
int GeometryRCB<Geometry,GlobalOrdinal>::getPointDestinationProc( 
    Teuchos::Array<double> coords ) const
{
    DTK_REQUIRE( 0 <= coords.size() && coords.size() <= 3 );
    DTK_REQUIRE( d_dimension == Teuchos::as<int>(coords.size()) );

    int proc = 0;
    DTK_REMEMBER( int zoltan_error );
#if HAVE_DTK_DBC
    zoltan_error = Zoltan_LB_Point_Assign( d_zz, &coords[0], &proc );
#else
    Zoltan_LB_Point_Assign( d_zz, &coords[0], &proc );
#endif
    DTK_CHECK( zoltan_error == ZOLTAN_OK );

    return proc;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the destination processes for a bounding box. This includes all
 * process domains that the box intersects.
 *
 * \param box The bounding box to get the destinations for.
 *
 * \return The GeometryRCB destination procs for the box
 */
template<class Geometry, class GlobalOrdinal>
Teuchos::Array<int>
GeometryRCB<Geometry,GlobalOrdinal>::getBoxDestinationProcs( 
    const BoundingBox& box ) const
{
    Teuchos::Tuple<double,6> box_bounds = box.getBounds();

    int num_procs = 0;
    Teuchos::Array<int> procs( d_comm->getSize() );

    DTK_REMEMBER( int zoltan_error );
#if HAVE_DTK_DBC
    zoltan_error = Zoltan_LB_Box_Assign( d_zz, 
					 box_bounds[0], box_bounds[1], 
					 box_bounds[2], box_bounds[3], 
					 box_bounds[4], box_bounds[5], 
					 &procs[0], &num_procs );
#else
    Zoltan_LB_Box_Assign( d_zz, 
			  box_bounds[0], box_bounds[1], box_bounds[2],
			  box_bounds[3], box_bounds[4], box_bounds[5], 
			  &procs[0], &num_procs );
#endif
    DTK_CHECK( zoltan_error == ZOLTAN_OK );

    procs.resize( num_procs );

    return procs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the number of centroids.
 */
template<class Geometry, class GlobalOrdinal>
int GeometryRCB<Geometry,GlobalOrdinal>::getNumberOfObjects( 
    void *data, int *ierr )
{
    RCP_GeometryManager geometry_manager = 
	*static_cast<RCP_GeometryManager*>( data );
    int num_geometry = 0;

    // We'll only count geometry if the geometry manager is not null.
    if ( !geometry_manager.is_null() )
    {
	Teuchos::ArrayView<short int> active_geom = 
	    geometry_manager->getActiveGeometry();

	num_geometry = 
	    std::accumulate( active_geom.begin(), active_geom.end(), 0 );
    }

    *ierr = ZOLTAN_OK;
    return num_geometry;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the local and global geometry ID's.
 */
template<class Geometry, class GlobalOrdinal>
void GeometryRCB<Geometry,GlobalOrdinal>::getObjectList( 
    void *data, int /*sizeGID*/, int /*sizeLID*/,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int /*wgt_dim*/, float * /*obj_wgts*/, int *ierr )
{
    RCP_GeometryManager geometry_manager = 
	*static_cast<RCP_GeometryManager*>( data );

    // We'll only build the geometry list is the geometry manager is not null.
    if ( !geometry_manager.is_null() )
    {
	// Note here that the local ID is being set to the geometry array
	// index.
	zoltan_id_type i = 0;
	GlobalOrdinal lid = 0;
	Teuchos::ArrayView<short int> active_geom =
		     geometry_manager->getActiveGeometry();
	Teuchos::ArrayView<short int>::const_iterator active_it;
	Teuchos::ArrayRCP<GlobalOrdinal> local_gids = 
	    geometry_manager->gids();
	typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator gid_it;
	typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator gids_begin = 
	    local_gids.begin();
	for ( gid_it = local_gids.begin(), active_it = active_geom.begin();
	      gid_it != local_gids.end(); 
	      ++gid_it, ++active_it )
	{
	    lid = std::distance( gids_begin, gid_it );

	    if ( *active_it )
	    {
		globalID[i] = static_cast<zoltan_id_type>( *gid_it );
		localID[i] = static_cast<zoltan_id_type>( lid );
		++i;
	    }
	}
    }

    *ierr = ZOLTAN_OK;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the dimension of the vertices.
 */
template<class Geometry, class GlobalOrdinal>
int GeometryRCB<Geometry,GlobalOrdinal>::getNumGeometry( void *data, int *ierr )
{
    int dimension = *static_cast<int*>( data );
    *ierr = ZOLTAN_OK;
    return dimension;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Zoltan callback for getting the centroid coordinates.
 */
template<class Geometry, class GlobalOrdinal>
void GeometryRCB<Geometry,GlobalOrdinal>::getGeometryList(
    void *data, int sizeGID, int sizeLID,
    int num_obj,
    ZOLTAN_ID_PTR /*globalID*/, ZOLTAN_ID_PTR /*localID*/,
    int num_dim, double *geom_vec, int *ierr )
{
    RCP_GeometryManager geometry_manager = 
	*static_cast<RCP_GeometryManager*>( data );

    // We will only supply centroid coordinates when the geometry exists.
    if ( !geometry_manager.is_null() )
    {
	Teuchos::ArrayRCP<Geometry> local_geometry = 
	    geometry_manager->geometry();
	typename Teuchos::ArrayRCP<Geometry>::const_iterator geom_it;
	Teuchos::ArrayView<short int> active_geom =
		     geometry_manager->getActiveGeometry();
	Teuchos::ArrayView<short int>::const_iterator active_it;

	// Check Zoltan for consistency.
	int geom_dim = geometry_manager->dim();
	DTK_CHECK( sizeGID == 1 );
	DTK_CHECK( sizeLID == 1 );
	DTK_CHECK( num_dim == Teuchos::as<int>(geom_dim) );

	if ( sizeGID != 1 || sizeLID != 1 || 
	     num_dim != Teuchos::as<int>(geom_dim) )
	{
	    *ierr = ZOLTAN_FATAL;
	    return;
	}
    
	// Zoltan needs interleaved coordinates.
	int n = 0;
	Teuchos::Array<double> geometry_coords;
	for ( geom_it = local_geometry.begin(), active_it = active_geom.begin();
	      geom_it != local_geometry.end();
	      ++geom_it, ++active_it )
	{
	    if ( *active_it )
	    {
		geometry_coords = GT::centroid( *geom_it );
		for ( int d = 0; d < geom_dim; ++d )
		{
		    geom_vec[ geom_dim*n + d ] = geometry_coords[d];
		}
		++n;
	    }
	}

	DTK_ENSURE( num_obj == n );
    }
	  
    *ierr = ZOLTAN_OK;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_GEOMETRYRCB_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_GeometryRCB_def.hpp
//---------------------------------------------------------------------------//
