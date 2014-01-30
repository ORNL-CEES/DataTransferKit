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
 * \file DTK_GeometryRendezvous_def.hpp
 * \author Stuart R. Slattery
 * \brief GeometryRendezvous definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_GEOMETRYRENDEZVOUS_DEF_HPP
#define DTK_GEOMETRYRENDEZVOUS_DEF_HPP

#include <set>
#include <algorithm>
#include <limits>

#include "DTK_DBC.hpp"
#include "DTK_CommIndexer.hpp"
#include "DTK_PartitionerFactory.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_as.hpp>

#include <Tpetra_Distributor.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param comm The communicator over which to build the rendezvous
 * decomposition. 
 *
 * \param dimension The dimension of the rendezvous decomposition. We need
 * this here because the geometry we get may or may not exist on every
 * process. This prevents global communications.
 *
 * \param global_box The global bounding box inside of which the rendezvous
 * decomposition will be generated.
 */
template<class Geometry, class GlobalOrdinal>
GeometryRendezvous<Geometry,GlobalOrdinal>::GeometryRendezvous( 
    const RCP_Comm& comm,
    const int dimension,
    const BoundingBox& global_box )
    : d_comm( comm )
    , d_dimension( dimension )
    , d_global_box( global_box )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Geometry, class GlobalOrdinal>
GeometryRendezvous<Geometry,GlobalOrdinal>::~GeometryRendezvous()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Build the rendezvous decomposition.
 *
 * \param geometry_manager The geometry to repartition to the rendezvous
 * decomposition. A null argument is valid here as a geometry may not exist in its
 * original decomposition on every process in the global communicator. It
 * will, however, be redistributed across every process in the global
 * communicator. 
 */
template<class Geometry, class GlobalOrdinal> 
void GeometryRendezvous<Geometry,GlobalOrdinal>::build( 
    const RCP_GeometryManager& geometry_manager )
{
    // Extract the geometry objects that are in the bounding box. These are
    // the pieces of the geometry that will be repartitioned.
    if ( !geometry_manager.is_null() ) 
    {
	getGeometryInBox( geometry_manager );
    }
    d_comm->barrier();

    // Construct the rendezvous partitioning for the geometry using the
    // vertices that are in the box.
    d_partitioner = PartitionerFactory::createGeometryPartitioner( 
	d_comm, geometry_manager, d_dimension );
    DTK_ENSURE( !d_partitioner.is_null() );
    d_partitioner->partition();

    // Send the geometry in the box to the rendezvous decomposition.
    sendGeometryToRendezvous( geometry_manager );
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Get the rendezvous destination processes for a blocked list of
 * coordinates that are in the primary decomposition.
 *
 * \param coords A blocked list of coordinates for which rendezvous
 * decomposition destination procs are desired.
 *
 * \return An array of the rendezvous decomposition destination procs. A proc
 * will be returned for each point in the same order as the points were
 * provided. 
 */
template<class Geometry, class GlobalOrdinal>
Teuchos::Array<int> 
GeometryRendezvous<Geometry,GlobalOrdinal>::procsContainingPoints(
    const Teuchos::ArrayRCP<double>& coords ) const
{
    Teuchos::Array<double> point( d_dimension );
    GlobalOrdinal num_points = coords.size() / d_dimension;
    Teuchos::Array<int> destination_procs( num_points );
    for ( GlobalOrdinal n = 0; n < num_points; ++n )
    {
	for ( int d = 0; d < d_dimension; ++d )
	{
	    point[d] = coords[ d*num_points + n ];
	}
	destination_procs[n] = d_partitioner->getPointDestinationProc( point );
    }

    return destination_procs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the rendezvous destination processes for a set of bounding
 * boxes. 
 *
 * \param boxes A list of boxes for which the rendezvous decomposition
 * destination procs are desired.
 *
 * \return A list of rendezvous decomposition procs for each box. A box may
 * have multiple procs that it spans.
 */
template<class Geometry, class GlobalOrdinal>
Teuchos::Array<Teuchos::Array<int> > 
GeometryRendezvous<Geometry,GlobalOrdinal>::procsContainingBoxes( 
    const Teuchos::Array<BoundingBox>& boxes ) const
{
    Teuchos::Array<Teuchos::Array<int> > box_procs( boxes.size() );
    Teuchos::Array<Teuchos::Array<int> >::iterator proc_iterator;
    Teuchos::Array<BoundingBox>::const_iterator box_iterator;
    for ( box_iterator = boxes.begin(), proc_iterator = box_procs.begin();
	  box_iterator != boxes.end();
	  ++box_iterator, ++proc_iterator )
    {
	*proc_iterator = d_partitioner->getBoxDestinationProcs( *box_iterator );
    }

    return box_procs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the native geometry gids in the rendezvous decomposition
 * containing a blocked list of coordinates also in the rendezvous
 * decomposition.
 * 
 * \param coords A blocked list of coordinates to search the geometry with. 

 * \param gids An array of the geometry ordinals the points were found
 * in. An element will be returned for each point in the order they were
 * provided in. If a point is not found in an geometry, return an invalid
 * geometry ordinal, std::numeric_limits<GlobalOrdinal>::max(), for that
 * point.
 *
 * \param geometry_src_procs The source procs that own the geometries. Once
 * proc is provided for each geometry in the order that the geometries were
 * provided. If a point is not found in an geometry, return an invalid
 * geometry source proc, -1, for that point.
 */
template<class Geometry, class GlobalOrdinal>
void GeometryRendezvous<Geometry,GlobalOrdinal>::geometryContainingPoints( 
    const Teuchos::ArrayRCP<double>& coords,
    Teuchos::Array<GlobalOrdinal>& gids,
    Teuchos::Array<int>& geometry_src_procs,
    const double geometric_tolerance ) const
{
    Teuchos::Array<double> point( d_dimension );
    GlobalOrdinal num_points = coords.size() / d_dimension;
    bool found_point = false;
    gids.resize( num_points );
    geometry_src_procs.resize( num_points );
    typename Teuchos::Array<Geometry>::const_iterator geom_it;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator gid_it;
    for ( GlobalOrdinal n = 0; n < num_points; ++n )
    {
	found_point = false;

	for ( int d = 0; d < d_dimension; ++d )
	{
	    point[d] = coords[ d*num_points + n ];
	}

	// Check for point inclusion in the geometry.
	for ( geom_it = d_rendezvous_geometry.begin(),
	       gid_it = d_rendezvous_gids.begin();
	      geom_it != d_rendezvous_geometry.end();
	      ++geom_it, ++gid_it )
	{
	    if ( !found_point )
	    {
		if ( GT::pointInGeometry(*geom_it, point, geometric_tolerance) )
		{
		    found_point = true;
		    gids[n] = *gid_it;
		    geometry_src_procs[n] = 
			d_geometry_src_procs_map.find( *gid_it )->second;
		}
	    }
	}

	// If we didnt find the point return an invalid geometry gid and
	// source proc.
	if ( !found_point )
	{
	    gids[n] = std::numeric_limits<GlobalOrdinal>::max();
	    geometry_src_procs[n] = -1;
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief For a list of geometry in the rendezvous decomposition, get their
 * source procs.
 *
 * \param gids The geometry gids to get source procs for.
 *
 * \return The source procs for the geometry.
 */
template<class Geometry, class GlobalOrdinal>
Teuchos::Array<int> 
GeometryRendezvous<Geometry,GlobalOrdinal>::geometrySourceProcs(
    const Teuchos::Array<GlobalOrdinal>& gids )
{
    Teuchos::Array<int> source_procs( gids.size() );
    Teuchos::Array<int>::iterator source_proc_iterator;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator gids_iterator;
    for ( gids_iterator = gids.begin(),
      source_proc_iterator = source_procs.begin();
	  gids_iterator != gids.end();
	  ++gids_iterator, ++source_proc_iterator )
    {
	*source_proc_iterator = 
	    d_geometry_src_procs_map.find( *gids_iterator )->second;
    }

    return source_procs;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Extract the geometry vertices and elements that are in a bounding
 * box. 
 *
 * \param geometry_manager The geometry to search the box with.
 */
template<class Geometry, class GlobalOrdinal>
void GeometryRendezvous<Geometry,GlobalOrdinal>::getGeometryInBox( 
    const RCP_GeometryManager& geometry_manager )
{
    // If the bounding box of a geometry intersects the rendezvous bounding
    // box then it is made active.
    Teuchos::ArrayRCP<Geometry> geometry = geometry_manager->geometry();
    typename Teuchos::ArrayRCP<Geometry>::const_iterator geometry_it;
    Teuchos::Array<short int> active_geometry( geometry.size(), 0 );
    Teuchos::Array<short int>::iterator active_it;
    BoundingBox dummy_box;

    for ( geometry_it = geometry.begin(), active_it = active_geometry.begin();
	  geometry_it != geometry.end();
	  ++geometry_it, ++active_it )
    {
	if ( BoundingBox::intersectBoxes( GT::boundingBox( *geometry_it ),
					  d_global_box, dummy_box ) )
	{
	    *active_it = 1;
	}
    }

    geometry_manager->setActiveGeometry( active_geometry );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Send the geometry to the rendezvous decomposition.
 *
 * \param geometry_manager The geometry to send to the rendezvous
 * decomposition. 
 */
template<class Geometry, class GlobalOrdinal>
void GeometryRendezvous<Geometry,GlobalOrdinal>::sendGeometryToRendezvous( 
    const RCP_GeometryManager& geometry_manager )
{
    // Set a value for geometry existence.
    bool geometry_exists = true;
    if ( geometry_manager.is_null() ) geometry_exists = false;

    // Setup a geometry data.
    Teuchos::ArrayRCP<Geometry> geometry(0);
    Teuchos::ArrayRCP<GlobalOrdinal> gids(0);
    RCP_Comm geometry_comm;
    Teuchos::Array<BoundingBox> geom_boxes(0);
    Teuchos::ArrayView<short int> active_geom;
    if ( geometry_exists )
    {
	geometry = geometry_manager->geometry();
	gids = geometry_manager->gids();
	geometry_comm = geometry_manager->comm();
	geom_boxes = geometry_manager->boundingBoxes();
	active_geom = geometry_manager->getActiveGeometry();
    }
    d_comm->barrier();
    CommIndexer geometry_indexer( d_comm, geometry_comm );

    // Get the rendezvous destination procs for the geometry bounding boxes.
    Teuchos::Array<Teuchos::Array<int> > box_procs =
	procsContainingBoxes( geom_boxes );

    // Unroll the geometry and procs into vectors for those that are active.
    Teuchos::Array<Teuchos::Array<int> >::const_iterator outer_it;
    Teuchos::Array<int>::const_iterator proc_it;
    typename Teuchos::ArrayRCP<Geometry>::const_iterator geom_it;
    typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator gid_it;
    Teuchos::ArrayView<short int>::const_iterator active_it;
    Teuchos::Array<Geometry> geom_to_rendezvous;
    Teuchos::Array<GlobalOrdinal> gids_to_rendezvous;
    Teuchos::Array<int> geom_rendezvous_procs;

    for( outer_it = box_procs.begin(), geom_it = geometry.begin(),
	   gid_it = gids.begin(), active_it = active_geom.begin();
	 outer_it != box_procs.end();
	 ++outer_it, ++geom_it, ++gid_it, ++active_it )
    {
	if ( *active_it )
	{
	    for ( proc_it = outer_it->begin();
		  proc_it != outer_it->end();
		  ++proc_it )
	    {
		geom_to_rendezvous.push_back( *geom_it );
		gids_to_rendezvous.push_back( *gid_it );
		geom_rendezvous_procs.push_back( *proc_it );
	    }
	}
    }

    // Distribute the geometry to the rendezvous decomposition.
    Tpetra::Distributor geometry_distributor( d_comm );
    GlobalOrdinal num_import_geom = 
	geometry_distributor.createFromSends( geom_rendezvous_procs() );
    geom_rendezvous_procs.clear();

    Teuchos::Array<Geometry> import_geom( num_import_geom );
    Teuchos::ArrayView<const Geometry> src_geom_view = geom_to_rendezvous();
    geometry_distributor.doPostsAndWaits( src_geom_view, 1, import_geom() );
    geom_to_rendezvous.clear();

    Teuchos::Array<GlobalOrdinal> import_gids( num_import_geom );
    Teuchos::ArrayView<const GlobalOrdinal> src_gids_view = 
	gids_to_rendezvous();
    geometry_distributor.doPostsAndWaits( src_gids_view, 1, import_gids() );
    gids_to_rendezvous.clear();

    // Extract the geometry source procs from the distributor.
    Teuchos::ArrayView<const int> from_images = 
	geometry_distributor.getImagesFrom();
    Teuchos::ArrayView<const std::size_t> from_lengths = 
	geometry_distributor.getLengthsFrom();
    Teuchos::Array<int> geometry_src_procs;
    for ( int i = 0; i < (int) from_images.size(); ++i )
    {
	for ( std::size_t j = 0; j < from_lengths[i]; ++j )
	{
	    geometry_src_procs.push_back( from_images[i] );
	}
    }
    DTK_CHECK( Teuchos::as<GlobalOrdinal>(geometry_src_procs.size()) 
		   == num_import_geom );

    // Build a unique set of local geometry, gids, and the source procs map.
    d_rendezvous_geometry.clear();
    d_rendezvous_gids.clear();
    std::set<GlobalOrdinal> import_gids_set;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator import_gids_it;
    typename Teuchos::Array<Geometry>::const_iterator import_geom_it;
    Teuchos::Array<int>::const_iterator geometry_src_procs_it;
    for ( import_gids_it = import_gids.begin(),
	  import_geom_it = import_geom.begin(),
       geometry_src_procs_it = geometry_src_procs.begin();
	  import_gids_it != import_gids.end();
	  ++import_gids_it, ++import_geom_it, ++geometry_src_procs_it )
    {
	if ( import_gids_set.insert( *import_gids_it ).second )
	{
	    d_geometry_src_procs_map[ *import_gids_it ] = 
		*geometry_src_procs_it;

	    d_rendezvous_geometry.push_back( *import_geom_it );

	    d_rendezvous_gids.push_back( *import_gids_it );
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_GEOMETRYRENDEZVOUS_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_GeometryRendezvous_def.hpp
//---------------------------------------------------------------------------//
