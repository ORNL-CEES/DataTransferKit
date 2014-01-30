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
 * \file DTK_VolumeSourceMap_def.hpp
 * \author Stuart R. Slattery
 * \brief Volume source map definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_VOLUMESOURCEMAP_DEF_HPP
#define DTK_VOLUMESOURCEMAP_DEF_HPP

#include <algorithm>
#include <limits>
#include <set>

#include "DTK_FieldTools.hpp"
#include "DTK_FieldTraits.hpp"
#include "DTK_DBC.hpp"
#include "DTK_GeometryRendezvous.hpp"
#include "DTK_BoundingBox.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_Import.hpp>
#include <Tpetra_Distributor.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Export.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param comm The communicator over which the map is generated.
 *
 * \param dimension The dimension of the map. This should be consistent with
 * all source and target objects (i.e. only 3 dimensional geometries will be
 * accepted with a 3 dimensional map). We need this here so we have a global
 * baseline for all objects that may or may not exist on all processes.
 * 
 * \param store_missed_points Set to true if it is desired to keep track of
 * the local target points missed during map generation. The default value is
 * false. 
 *
 * \param geometric_tolerance Tolerance used for point-in-geometry checks. The
 * default value is 1.0e-6.
 */
template<class Geometry, class GlobalOrdinal, class CoordinateField>
VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::VolumeSourceMap(
    const RCP_Comm& comm, const int dimension, 
    bool store_missed_points,
    const double geometric_tolerance )
    : d_comm( comm )
    , d_dimension( dimension )
    , d_store_missed_points( store_missed_points )
    , d_geometric_tolerance( geometric_tolerance )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Geometry, class GlobalOrdinal, class CoordinateField>
VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::~VolumeSourceMap()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Generate the volume source map.
 *
 * \param source_geometry_manager Source geometry in the volume source
 * problem. A null RCP is a valid argument. This will be the case when a
 * geometry manager is only constructed on a subset of the processes that the
 * shared domain map is constructed over. Note that the source geometry must
 * exist only on processes that reside within the VolumeSourceMap
 * communicator.
 *
 * \param target_coord_manager Target coordinates in the shared domain
 * problem. A null RCP is a valid argument. This will be the case when a field
 * manager is only constructed on a subset of the processes that the shared
 * domain map is constructed over. Note that the target coordinates must exist
 * only on processes that reside within the VolumeSourceMap communicator.
 */
template<class Geometry, class GlobalOrdinal, class CoordinateField>
void VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::setup( 
    const RCP_GeometryManager& source_geometry_manager, 
    const RCP_CoordFieldManager& target_coord_manager )
{
    // Create existence values for the managers.
    bool source_exists = true;
    if ( source_geometry_manager.is_null() ) source_exists = false;
    bool target_exists = true;
    if ( target_coord_manager.is_null() ) target_exists = false;

    // Create local to global process indexers for the managers.
    RCP_Comm source_comm;
    if ( source_exists )
    {
	source_comm = source_geometry_manager->comm();
    }
    d_source_indexer = CommIndexer( d_comm, source_comm );
    RCP_Comm target_comm;
    if ( target_exists )
    {
	target_comm = target_coord_manager->comm();
    }
    d_target_indexer = CommIndexer( d_comm, target_comm );

    // Check the source and target dimensions for consistency.
    if ( source_exists )
    {
	DTK_REQUIRE( source_geometry_manager->dim() == d_dimension );
    }
    if ( target_exists )
    {
	DTK_REQUIRE( CFT::dim( *target_coord_manager->field() ) 
			  == d_dimension );
    }

    // Compute a unique global ordinal for each point in the coordinate field.
    Teuchos::Array<GlobalOrdinal> target_ordinals;
    computePointOrdinals( target_coord_manager, target_ordinals );

    // Build the data import map from the point global ordinals.
    Teuchos::ArrayView<const GlobalOrdinal> import_ordinal_view =
	target_ordinals();
    d_target_map = Tpetra::createNonContigMap<int,GlobalOrdinal>(
	import_ordinal_view, d_comm );
    DTK_ENSURE( !d_target_map.is_null() );

    // Get the global bounding box for the geometry.
    BoundingBox source_box;
    if ( source_exists )
    {
	source_box = source_geometry_manager->globalBoundingBox();
    }
    Teuchos::broadcast<int,BoundingBox>(
	*d_comm, d_source_indexer.l2g(0), 
	Teuchos::Ptr<BoundingBox>(&source_box) );

    // Get the global bounding box for the coordinate field.
    BoundingBox target_box;
    if ( target_exists )
    {
	target_box = FieldTools<CoordinateField>::coordGlobalBoundingBox(
	    *target_coord_manager->field(), target_coord_manager->comm() );
    }
    Teuchos::broadcast<int,BoundingBox>( 
	*d_comm, d_target_indexer.l2g(0), 
	Teuchos::Ptr<BoundingBox>(&target_box) );

    // Intersect the boxes to get the shared domain bounding box.
    BoundingBox shared_domain_box;
    bool has_intersect = BoundingBox::intersectBoxes( source_box, target_box, 
						      shared_domain_box );
    DTK_INSIST( has_intersect );

    // Build a rendezvous decomposition with the source geometry.
    GeometryRendezvous<Geometry,GlobalOrdinal> rendezvous( 
	d_comm, d_dimension, shared_domain_box );
    rendezvous.build( source_geometry_manager );

    // Determine the rendezvous destination proc of each point in the
    // coordinate field.
    Teuchos::ArrayRCP<double> coords_view(0,0.0);
    int coord_dim;
    if ( target_exists )
    {
	coord_dim = CFT::dim( *target_coord_manager->field() );
	coords_view = FieldTools<CoordinateField>::nonConstView( 
	    *target_coord_manager->field() );
    }
    Teuchos::broadcast<int,int>( 
	*d_comm, d_target_indexer.l2g(0), Teuchos::Ptr<int>(&coord_dim) );

    Teuchos::Array<int> rendezvous_procs = 
	rendezvous.procsContainingPoints( coords_view );

    // Get the target points that are in the box in which the rendezvous
    // decomposition was generated. 
    Teuchos::Array<GlobalOrdinal> targets_in_box;
    if ( target_exists )
    {
	getTargetPointsInBox( rendezvous.getBox(), 
			      *target_coord_manager->field(),
			      target_ordinals, targets_in_box );
    }

    // Extract those target points that are not in the box. We don't want to
    // send these to the rendezvous decomposition.
    Teuchos::Array<GlobalOrdinal> not_in_box;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator in_box_iterator;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator in_box_begin =
	targets_in_box.begin();
    for ( in_box_iterator = targets_in_box.begin();
	  in_box_iterator != targets_in_box.end();
	  ++in_box_iterator )
    {
	if ( *in_box_iterator == std::numeric_limits<GlobalOrdinal>::max() )
	{
	    not_in_box.push_back( 
		std::distance( in_box_begin, in_box_iterator ) );
	}
    }
    std::reverse( not_in_box.begin(), not_in_box.end() );

    typename Teuchos::Array<GlobalOrdinal>::const_iterator not_in_box_iterator;
    for ( not_in_box_iterator = not_in_box.begin();
	  not_in_box_iterator != not_in_box.end();
	  ++not_in_box_iterator )
    {
	rendezvous_procs.remove( *not_in_box_iterator );
    }
    not_in_box.clear();

    typename Teuchos::Array<GlobalOrdinal>::iterator targets_bound =
	std::remove( targets_in_box.begin(), targets_in_box.end(), 
		     std::numeric_limits<GlobalOrdinal>::max() );
    GlobalOrdinal targets_in_box_size = 
	std::distance( targets_in_box.begin(), targets_bound );

    targets_in_box.resize( targets_in_box_size );
    rendezvous_procs.resize( targets_in_box_size );

    // Via an inverse communication operation, move the global point ordinals
    // that are in the rendezvous decomposition box to the rendezvous
    // decomposition.
    Teuchos::ArrayView<const GlobalOrdinal> targets_in_box_view = 
	targets_in_box();
    Tpetra::Distributor target_to_rendezvous_distributor( d_comm );
    GlobalOrdinal num_rendezvous_points = 
	target_to_rendezvous_distributor.createFromSends( rendezvous_procs() );
    Teuchos::Array<GlobalOrdinal> rendezvous_points( num_rendezvous_points );
    target_to_rendezvous_distributor.doPostsAndWaits( 
	targets_in_box_view, 1, rendezvous_points() );

    // Setup target-to-rendezvous communication.
    Teuchos::ArrayView<const GlobalOrdinal> rendezvous_points_view =
	rendezvous_points();
    RCP_TpetraMap rendezvous_coords_map = 
      Tpetra::createNonContigMap<int,GlobalOrdinal>(
        rendezvous_points_view, d_comm );
    Tpetra::Export<int,GlobalOrdinal> 
      target_to_rendezvous_exporter( d_target_map, rendezvous_coords_map );

    // Move the target coordinates to the rendezvous decomposition.
    GlobalOrdinal num_points = target_ordinals.size();
    Teuchos::RCP< Tpetra::MultiVector<double,int,GlobalOrdinal> > 
	target_coords =	Tpetra::createMultiVectorFromView( 
	    d_target_map, coords_view, num_points, coord_dim );
    Tpetra::MultiVector<double,int,GlobalOrdinal> rendezvous_coords( 
	rendezvous_coords_map, coord_dim );
    rendezvous_coords.doExport( *target_coords, target_to_rendezvous_exporter, 
				Tpetra::INSERT );

    // Search the rendezvous decomposition with the target points to get the
    // source geometry that contains them.
    Teuchos::Array<GlobalOrdinal> rendezvous_geometry;
    Teuchos::Array<int> rendezvous_geometry_src_procs;
    rendezvous.geometryContainingPoints( rendezvous_coords.get1dViewNonConst(),
					 rendezvous_geometry,
					 rendezvous_geometry_src_procs,
					 d_geometric_tolerance );

    // Get the points that were not in the geometry. If we're keeping track of
    // missed points, also make a list of those ordinals.
    GlobalOrdinal target_index;
    Teuchos::Array<GlobalOrdinal> not_in_geometry;
    Teuchos::Array<GlobalOrdinal> missed_in_geometry_idx;
    Teuchos::Array<GlobalOrdinal> missed_in_geometry_ordinal;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
	rendezvous_geometry_iterator;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
	rendezvous_geometry_begin = rendezvous_geometry.begin();
    for ( rendezvous_geometry_iterator = rendezvous_geometry.begin();
	  rendezvous_geometry_iterator != rendezvous_geometry.end();
	  ++rendezvous_geometry_iterator )
    {
	if ( *rendezvous_geometry_iterator == 
	     std::numeric_limits<GlobalOrdinal>::max() )
	{
	    target_index = std::distance( rendezvous_geometry_begin, 
					  rendezvous_geometry_iterator );

	    not_in_geometry.push_back( target_index );

	    if ( d_store_missed_points )
	    {
		missed_in_geometry_idx.push_back( target_index );
		missed_in_geometry_ordinal.push_back( 
		    rendezvous_points[target_index] );
	    }
	}
    }

    // If we're keeping track of missed points, send their global ordinals
    // back to the target decomposition so that we can add them to the list.
    if ( d_store_missed_points )
    {
	// Extract the missed point target procs from the target-to-rendezvous
	// distributor.
	Teuchos::ArrayView<const int> from_images = 
	    target_to_rendezvous_distributor.getImagesFrom();
	Teuchos::ArrayView<const std::size_t> from_lengths = 
	    target_to_rendezvous_distributor.getLengthsFrom();
	Teuchos::Array<int> point_target_procs;
	for ( int i = 0; i < (int) from_images.size(); ++i )
	{
	    for ( std::size_t j = 0; j < from_lengths[i]; ++j )
	    {
		point_target_procs.push_back( from_images[i] );
	    }
	}
	DTK_CHECK( Teuchos::as<GlobalOrdinal>(point_target_procs.size()) ==
		       num_rendezvous_points );

	// Build a list of target procs for the missed points.
	Teuchos::Array<int> missed_target_procs( missed_in_geometry_idx.size() );
	for ( int n = 0; n < (int) missed_in_geometry_idx.size(); ++n )
	{
	    missed_target_procs[n] = 
		point_target_procs[ missed_in_geometry_idx[n] ];
	}
	point_target_procs.clear();

	// Send the missed points back to the target decomposition through an
	// inverse communication operation and add them to the list.
	Teuchos::ArrayView<const GlobalOrdinal> missed_in_geometry_ordinal_view = 
	    missed_in_geometry_ordinal();
	Tpetra::Distributor target_to_rendezvous_distributor( d_comm );
	GlobalOrdinal num_missed_targets = 
	    target_to_rendezvous_distributor.createFromSends( 
		missed_target_procs() );
	GlobalOrdinal offset = d_missed_points.size();
	d_missed_points.resize( offset + num_missed_targets );
	target_to_rendezvous_distributor.doPostsAndWaits( 
	    missed_in_geometry_ordinal_view, 1, 
	    d_missed_points.view( offset, num_missed_targets ) );

	// Convert the missed point global indices to local indices and add
	// them to the list.
	for ( GlobalOrdinal n = offset; n < offset+num_missed_targets; ++n )
	{
	    d_missed_points[n] = 
		d_target_g2l.find( d_missed_points[n] )->second;
	}
    }
    missed_in_geometry_idx.clear();
    missed_in_geometry_ordinal.clear();

    // Extract the points we didn't find in any geometry in the rendezvous
    // decomposition and their corresponding geometry. We don't want to send
    // these to the source.
    std::reverse( not_in_geometry.begin(), not_in_geometry.end() );
    typename Teuchos::Array<GlobalOrdinal>::const_iterator not_in_geometry_iterator;
    for ( not_in_geometry_iterator = not_in_geometry.begin();
	  not_in_geometry_iterator != not_in_geometry.end();
	  ++not_in_geometry_iterator )
    {
	rendezvous_points.remove( *not_in_geometry_iterator );
    }
    not_in_geometry.clear();

    typename Teuchos::Array<GlobalOrdinal>::iterator rendezvous_geometry_bound =
	std::remove( rendezvous_geometry.begin(), 
		     rendezvous_geometry.end(), 
		     std::numeric_limits<GlobalOrdinal>::max() );
    GlobalOrdinal rendezvous_geometry_size = 
	std::distance( rendezvous_geometry.begin(), rendezvous_geometry_bound );

    typename Teuchos::Array<int>::iterator rendezvous_geometry_src_procs_bound =
	std::remove( rendezvous_geometry_src_procs.begin(), 
		     rendezvous_geometry_src_procs.end(), -1 );
    GlobalOrdinal rendezvous_geometry_src_procs_size = 
	std::distance( rendezvous_geometry_src_procs.begin(), 
		       rendezvous_geometry_src_procs_bound );
    DTK_CHECK( rendezvous_geometry_size == 
		   rendezvous_geometry_src_procs_size );

    rendezvous_geometry.resize( rendezvous_geometry_size );
    rendezvous_geometry_src_procs.resize( rendezvous_geometry_src_procs_size );
    rendezvous_points.resize( rendezvous_geometry_size );

    // Setup rendezvous-to-source distributor.
    Tpetra::Distributor rendezvous_to_src_distributor( d_comm );
    GlobalOrdinal num_source_geometry = 
	rendezvous_to_src_distributor.createFromSends( 
	    rendezvous_geometry_src_procs() );

    // Send the rendezvous geometry to the source decomposition via inverse
    // communication. 
    Teuchos::ArrayView<const GlobalOrdinal> rendezvous_geometry_view =
	rendezvous_geometry();
    d_source_geometry.resize( num_source_geometry );
    rendezvous_to_src_distributor.doPostsAndWaits( rendezvous_geometry_view, 1,
						   d_source_geometry() );

    // Send the rendezvous point global ordinals to the source decomposition
    // via inverse communication.
    Teuchos::ArrayView<const GlobalOrdinal> reduced_rendezvous_points_view =
	rendezvous_points();
    Teuchos::Array<GlobalOrdinal> source_points( num_source_geometry );
    rendezvous_to_src_distributor.doPostsAndWaits( 
	reduced_rendezvous_points_view, 1, source_points() );

    // Build the source map from the target ordinals.
    Teuchos::ArrayView<const GlobalOrdinal> source_points_view = 
	source_points();
    d_source_map = Tpetra::createNonContigMap<int,GlobalOrdinal>( 
	source_points_view, d_comm );
    DTK_ENSURE( !d_source_map.is_null() );

    // Send the rendezvous point coordinates to the source decomposition.
    Tpetra::Export<int,GlobalOrdinal> rendezvous_to_source_exporter( 
	rendezvous_coords_map, d_source_map );
    d_target_coords.resize( num_source_geometry*coord_dim );
    Teuchos::RCP< Tpetra::MultiVector<double,int,GlobalOrdinal> >
	source_coords = Tpetra::createMultiVectorFromView( 
	    d_source_map, Teuchos::arcpFromArray( d_target_coords ), 
	    num_source_geometry, coord_dim );
    source_coords->doExport( rendezvous_coords, rendezvous_to_source_exporter,
			     Tpetra::INSERT );

    // Build the source-to-target importer.
    d_source_to_target_importer = 
      Teuchos::rcp( new Tpetra::Import<int,GlobalOrdinal>(
          d_source_map, d_target_map ) );
    DTK_ENSURE( !d_source_to_target_importer.is_null() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply the volume source map for a valid source field evaluator and
 * target data space to the target points that were mapped.
 *
 * \param source_evaluator Function evaluator used to apply the mapping. This
 * FieldEvaluator must be valid for the source geometry used to generate the
 * map.
 *
 * \param target_space_manager Target space into which the function
 * evaluations will be written. Enough space must be allocated to hold
 * evaluations at all points in all dimensions of the field.
 */
template<class Geometry, class GlobalOrdinal, class CoordinateField>
template<class SourceField, class TargetField>
void VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::apply( 
    const Teuchos::RCP< FieldEvaluator<GlobalOrdinal,SourceField> >& source_evaluator,
    Teuchos::RCP< FieldManager<TargetField> >& target_space_manager )
{
    typedef FieldTraits<SourceField> SFT;
    typedef FieldTraits<TargetField> TFT;

    // Set existence values for the source and target.
    bool source_exists = true;
    if ( source_evaluator.is_null() ) source_exists = false;
    bool target_exists = true;
    if ( target_space_manager.is_null() ) target_exists = false;

    // Evaluate the source function at the target points and construct a view
    // of the function evaluations.
    int source_dim;
    Teuchos::ArrayRCP<typename SFT::value_type> source_field_copy(0,0);
    if ( source_exists )
    {
	SourceField function_evaluations = 
	    source_evaluator->evaluate( 
		Teuchos::arcpFromArray( d_source_geometry ),
		Teuchos::arcpFromArray( d_target_coords ) );

	source_dim = SFT::dim( function_evaluations );

	source_field_copy =    
	    FieldTools<SourceField>::copy( function_evaluations );
    }
    Teuchos::broadcast<int,int>( *d_comm, d_source_indexer.l2g(0),
				 Teuchos::Ptr<int>(&source_dim) );

    // Build a multivector for the function evaluations.
    GlobalOrdinal source_size = source_field_copy.size() / source_dim;
    Teuchos::RCP<Tpetra::MultiVector<typename SFT::value_type, int, GlobalOrdinal> > 
	source_vector = Tpetra::createMultiVectorFromView( 
	    d_source_map, source_field_copy, source_size, source_dim );

    // Construct a view of the target space.
    int target_dim;
    Teuchos::ArrayRCP<typename TFT::value_type> target_field_view(0,0);
    if ( target_exists )
    {
	target_field_view = FieldTools<TargetField>::nonConstView( 
	    *target_space_manager->field() );

	target_dim = TFT::dim( *target_space_manager->field() );
    }
    Teuchos::broadcast<int,int>( *d_comm, d_target_indexer.l2g(0),
				 Teuchos::Ptr<int>(&target_dim) );
    
    // Check that the source and target have the same field dimension.
    DTK_REQUIRE( source_dim == target_dim );

    // Verify that the target space has the proper amount of memory allocated.
    GlobalOrdinal target_size = target_field_view.size() / target_dim;
    if ( target_exists )
    {
	DTK_REQUIRE( 
	    target_size == Teuchos::as<GlobalOrdinal>(
		d_target_map->getNodeNumElements()) );
    }
    
    // Build a multivector for the target space.
    Teuchos::RCP<Tpetra::MultiVector<typename TFT::value_type, int, GlobalOrdinal> > 
	target_vector =	Tpetra::createMultiVectorFromView( 
	    d_target_map, target_field_view, target_size, target_dim );

    // Fill the target space with zeros so that points we didn't map get some
    // data.
    if ( target_exists )
    {
	FieldTools<TargetField>::putScalar( 
	    *target_space_manager->field(), 0.0 );
    }

    // Move the data from the source decomposition to the target
    // decomposition.
    target_vector->doImport( *source_vector, *d_source_to_target_importer, 
			     Tpetra::INSERT );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the points missed in the map generation.
 *
 * \return If store_missed_points is true, return the local indices of the
 *  points provided by target_coord_manager that were not mapped. An exception
 *  will be thrown if store_missed_points is false. Returns a null view if all
 *  points have been mapped or the map has not yet been generated.
*/
template<class Geometry, class GlobalOrdinal, class CoordinateField>
Teuchos::ArrayView<const GlobalOrdinal>
VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::getMissedTargetPoints() const
{
    DTK_REQUIRE( d_store_missed_points );
    
    return d_missed_points();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the points missed in the map generation.
 *
 * \return If store_missed_points is true, return the local indices of the
 *  points provided by target_coord_manager that were not mapped. An exception
 *  will be thrown if store_missed_points is false. Returns a null view if all
 *  points have been mapped or the map has not yet been generated.
*/
template<class Geometry, class GlobalOrdinal, class CoordinateField>
Teuchos::ArrayView<GlobalOrdinal> 
VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::getMissedTargetPoints()
{
    DTK_REQUIRE( d_store_missed_points );
    
    return d_missed_points();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute globally unique ordinals for the target points. Here an
 * invalid ordinal will be designated as the maximum value as specified by the
 * limits header for the ordinal type. We do this so that 0 may be a valid
 * ordinal.
 *
 * \param target_coords The coordinates to compute global ordinals for.
 *
 * \param target_ordinals The computed globally unique ordinals for the target
 * coordinates. 
 */
template<class Geometry, class GlobalOrdinal, class CoordinateField>
void 
VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::computePointOrdinals(
    const RCP_CoordFieldManager& target_coord_manager,
    Teuchos::Array<GlobalOrdinal>& target_ordinals )
{
    // Set an existence value for the target coords.
    bool target_exists = true;
    if ( target_coord_manager.is_null() ) target_exists = false;
    int comm_rank = d_comm->getRank();
    GlobalOrdinal local_size = 0;

    if ( target_exists )
    {
	int point_dim = CFT::dim( *target_coord_manager->field() );
	local_size = std::distance( 
	    CFT::begin( *target_coord_manager->field() ),
	    CFT::end( *target_coord_manager->field() ) ) / point_dim;
    }

    GlobalOrdinal global_max;
    Teuchos::reduceAll<int,GlobalOrdinal>( *d_comm,
					   Teuchos::REDUCE_MAX,
					   1,
					   &local_size,
					   &global_max );

    target_ordinals.resize( local_size );
    for ( GlobalOrdinal n = 0; n < local_size; ++n )
    {
	target_ordinals[n] = comm_rank*global_max + n;

	// If we're keeping track of missed points, we also need to build the
	// global-to-local ordinal map.
	if ( d_store_missed_points )
	{
	    d_target_g2l[ target_ordinals[n] ] = n;
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the target points that are in the rendezvous decomposition box.
 *
 * \param target_coord_manager The manager containing the target coordinates
 * to search the box with.
 *
 * \param target_ordinals The globally unique ordinals for the target
 * coordinates. 
 *
 * \param targets_in_box The global ordinals of the target coordinates in the
 * box. If a target point was not found in the box, return an invalid ordinal,
 * std::numeric_limits<GlobalOrdinal>::max(), in its positition.
 */
template<class Geometry, class GlobalOrdinal, class CoordinateField>
void VolumeSourceMap<Geometry,GlobalOrdinal,CoordinateField>::getTargetPointsInBox(
    const BoundingBox& box, const CoordinateField& target_coords,
    const Teuchos::Array<GlobalOrdinal>& target_ordinals,
    Teuchos::Array<GlobalOrdinal>& targets_in_box )
{
    // Expand the box in all directions by the geometric tolerance. Doing this
    // catches a few corner cases.
    Teuchos::Tuple<double,6> box_bounds = box.getBounds();
    for ( int d = 0; d < d_dimension; ++d )
    {
	box_bounds[d] -= d_geometric_tolerance;
	box_bounds[d+3] += d_geometric_tolerance;
    }
    BoundingBox expanded_box( box_bounds );

    Teuchos::ArrayRCP<const double> target_coords_view =
	FieldTools<CoordinateField>::view( target_coords );
    GlobalOrdinal dim_size = 
	FieldTools<CoordinateField>::dimSize( target_coords );

    DTK_REQUIRE( dim_size == 
		      Teuchos::as<GlobalOrdinal>(target_ordinals.size()) );

    targets_in_box.resize( dim_size );
    int field_dim = CFT::dim( target_coords );
    Teuchos::Array<double> target_point( field_dim );
    for ( GlobalOrdinal n = 0; n < dim_size; ++n )
    {
	for ( int d = 0; d < field_dim; ++d )
	{
	    target_point[d] = target_coords_view[ dim_size*d + n ];
	}

	if ( expanded_box.pointInBox( target_point ) )
	{
	    targets_in_box[n] = target_ordinals[n];
	}
	else
	{
	    targets_in_box[n] = std::numeric_limits<GlobalOrdinal>::max();
	}

	// If we're keeping track of the points not being mapped, add this
	// point's local index to the list if its not in the box.
	if ( d_store_missed_points && targets_in_box[n] == 
	     std::numeric_limits<GlobalOrdinal>::max() )
	{
	    d_missed_points.push_back(n);
	}
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_VOLUMESOURCEMAP_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_VolumeSourceMap_def.hpp
//---------------------------------------------------------------------------//


