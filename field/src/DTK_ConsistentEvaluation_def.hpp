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
 * \file DTK_ConsistentEvaluation.hpp
 * \author Stuart R. Slattery
 * \brief Consistent evaluation mapping definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CONSISTENTEVALUATION_DEF_HPP
#define DTK_CONSISTENTEVALUATION_DEF_HPP

#include <algorithm>
#include <set>

#include "DTK_FieldTools.hpp"
#include <DTK_Exception.hpp>
#include <DTK_Rendezvous.hpp>
#include <DTK_MeshTools.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_Import.hpp>
#include <Tpetra_Distributor.hpp>
#include <Tpetra_MultiVector.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Mesh, class CoordinateField>
ConsistentEvaluation<Mesh,CoordinateField>::ConsistentEvaluation( 
    const RCP_Comm& comm )
    : d_comm( comm )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Mesh, class CoordinateField>
ConsistentEvaluation<Mesh,CoordinateField>::~ConsistentEvaluation()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Setup for evaluation.
 */
template<class Mesh, class CoordinateField>
void ConsistentEvaluation<Mesh,CoordinateField>::setup( 
    const RCP_MeshManager& mesh_manager, 
    const CoordinateField& coordinate_field )
{
    // Get the global bounding box for the mesh.
    BoundingBox mesh_box = mesh_manager->globalBoundingBox();

    // Get the global bounding box for the coordinate field.
    BoundingBox coord_box = FieldTools<CoordinateField>::coordGlobalBoundingBox(
	coordinate_field, d_comm );

    // Intersect the boxes to get the rendezvous bounding box.
    BoundingBox rendezvous_box;
    bool has_intersect = 
	BoundingBox::intersectBoxes( mesh_box, coord_box, rendezvous_box );
    if ( !has_intersect )
    {
	throw MeshException( 
	    "Mesh and coordinate field domains do not intersect." );
    }

    // Build a rendezvous decomposition with the source mesh.
    Rendezvous<Mesh> rendezvous( d_comm, rendezvous_box );
    rendezvous.build( mesh_manager );

    // Compute a unique global ordinal for each point in the coordinate field.
    Teuchos::Array<GlobalOrdinal> point_ordinals = 
	computePointOrdinals( coordinate_field );

    // Build the data import map from the point global ordinals.
    Teuchos::ArrayView<const GlobalOrdinal> import_ordinal_view =
	point_ordinals();
    d_import_map = Tpetra::createNonContigMap<GlobalOrdinal>(
	import_ordinal_view, d_comm );
    testPostcondition( d_import_map != Teuchos::null,
		       "Error creating data import map." );

    // Determine the rendezvous destination proc of each point in the
    // coordinate field.
    std::size_t coord_dim = CFT::dim( coordinate_field );
    typename CFT::size_type num_coords = CFT::size( coordinate_field );
    Teuchos::ArrayRCP<typename CFT::value_type> coords_view;
    if ( num_coords == 0 )
    {
	coords_view = Teuchos::ArrayRCP<typename CFT::value_type>( 0, 0.0 );
    }
    else
    {
	coords_view = 
	    FieldTools<CoordinateField>::nonConstView( coordinate_field );
    }
    Teuchos::Array<int> rendezvous_procs = 
	rendezvous.getRendezvousProcs( coords_view );

    // Via an inverse communication operation, move the global point ordinals
    // to the rendezvous decomposition.
    Tpetra::Distributor point_distributor( d_comm );
    GlobalOrdinal num_rendezvous_points = 
	point_distributor.createFromSends( rendezvous_procs() );
    Teuchos::Array<GlobalOrdinal> 
	rendezvous_points( num_rendezvous_points );
    point_distributor.doPostsAndWaits( 
	import_ordinal_view, 1, rendezvous_points() );

    // Setup Target-to-rendezvous communication.
    Teuchos::ArrayView<const GlobalOrdinal> rendezvous_points_view =
	rendezvous_points();
    RCP_TpetraMap rendezvous_import_map = 
	Tpetra::createNonContigMap<GlobalOrdinal>(
	    rendezvous_points_view, d_comm );
    Tpetra::Export<GlobalOrdinal> 
	point_exporter( d_import_map, rendezvous_import_map );

    // Move the target coordinates to the rendezvous decomposition.
    GlobalOrdinal num_points = num_coords / coord_dim;
    Teuchos::RCP< Tpetra::MultiVector<double,GlobalOrdinal> > 
	target_coords =	createMultiVectorFromView( d_import_map, coords_view,
						   num_points, coord_dim );
    Tpetra::MultiVector<double,GlobalOrdinal> rendezvous_coords( 
	rendezvous_import_map, coord_dim );
    rendezvous_coords.doExport( *target_coords, point_exporter, 
				Tpetra::INSERT );

    // Search the rendezvous decomposition with the target points to get the
    // source elements that contain them.
    Teuchos::Array<GlobalOrdinal> rendezvous_elements =
	rendezvous.getElements( rendezvous_coords.get1dViewNonConst() );

    // Build a unique list of rendezvous elements.
    Teuchos::Array<GlobalOrdinal> rendezvous_element_set =
	rendezvous_elements;
    std::sort( rendezvous_element_set.begin(), rendezvous_element_set.end() );
    typename Teuchos::Array<GlobalOrdinal>::iterator unique_bound =
	std::unique( rendezvous_element_set.begin(), 
		     rendezvous_element_set.end() );
    rendezvous_element_set.resize( 
	unique_bound - rendezvous_element_set.begin() );

    // Setup source-to-rendezvous communication.
    Teuchos::ArrayView<const GlobalOrdinal> rendezvous_elem_set_view =
	rendezvous_element_set();
    RCP_TpetraMap rendezvous_element_map = 
	Tpetra::createNonContigMap<GlobalOrdinal>( 
	    rendezvous_elem_set_view, d_comm );

    Teuchos::Array<GlobalOrdinal> 
	mesh_elements( mesh_manager->localNumElements() );
    BlockIterator block_iterator;
    for ( block_iterator = mesh_manager->blocksBegin();
	  block_iterator != mesh_manager->blocksEnd();
	  ++block_iterator )
    {
	Teuchos::ArrayRCP<const GlobalOrdinal> mesh_element_arcp =
	    MeshTools<Mesh>::elementsView( mesh );
    }

    Teuchos::ArrayView<const GlobalOrdinal> mesh_element_view = 
	mesh_elements();
    RCP_TpetraMap mesh_element_map = 
	Tpetra::createNonContigMap<GlobalOrdinal>( 
	    mesh_element_view, d_comm );
    mesh_elements.clear();

    Tpetra::Import<GlobalOrdinal> source_importer( 
	mesh_element_map, rendezvous_element_map );

    // Elements send their source decomposition proc to the rendezvous
    // decomposition.
    Tpetra::MultiVector<int,GlobalOrdinal> 
	source_procs( mesh_element_map, 1 );
    source_procs.putScalar( d_comm->getRank() );
    Tpetra::MultiVector<int,GlobalOrdinal> 
	rendezvous_source_procs( rendezvous_element_map, 1 );
    rendezvous_source_procs.doImport( source_procs, source_importer, 
				      Tpetra::INSERT );

    // Set the destination procs for the rendezvous-to-source communication.
    Teuchos::ArrayRCP<const int> rendezvous_source_procs_view =
	rendezvous_source_procs.get1dView();
    GlobalOrdinal num_rendezvous_elements = rendezvous_elements.size();
    GlobalOrdinal dest_proc_idx;
    Teuchos::Array<int> 
	target_destinations( num_rendezvous_elements );
    for ( int n = 0; n < num_rendezvous_elements; ++n )
    {
	dest_proc_idx = std::distance( 
	    rendezvous_element_set.begin(),
	    std::find( rendezvous_element_set.begin(),
		       rendezvous_element_set.end(),
		       rendezvous_elements[n] ) );
	target_destinations[n] = rendezvous_source_procs_view[ dest_proc_idx ];
    }
    rendezvous_element_set.clear();

    // Send the rendezvous elements to the source decomposition via inverse
    // communication. 
    Teuchos::ArrayView<const GlobalOrdinal> rendezvous_elements_view =
	rendezvous_elements();
    Tpetra::Distributor source_distributor( d_comm );
    GlobalOrdinal num_source_elements = 
	source_distributor.createFromSends( target_destinations() );
    d_source_elements.resize( num_source_elements );
    source_distributor.doPostsAndWaits( rendezvous_elements_view, 1,
					d_source_elements() );

    // Send the rendezvous point coordinates to the source decomposition via
    // inverse communication. 
    Teuchos::ArrayView<const double> rendezvous_coords_view =
	rendezvous_coords.get1dView()();
    d_target_coords.resize( num_source_elements*coord_dim );
    source_distributor.doPostsAndWaits( rendezvous_coords_view, coord_dim,
					d_target_coords() );

    // Send the rendezvous point global ordinals to the source decomposition
    // via inverse communication.
    Teuchos::Array<GlobalOrdinal> source_points( num_source_elements );
    source_distributor.doPostsAndWaits( rendezvous_points_view, 1,
					source_points() );

    // Build the data export map from the coordinate ordinals as well as
    // populate the source element/target coordinate pairings.
    Teuchos::ArrayView<const GlobalOrdinal> source_points_view =
	source_points();
    d_export_map = Tpetra::createNonContigMap<GlobalOrdinal>(
	source_points_view, d_comm );
    testPostcondition( d_export_map != Teuchos::null,
		       "Error creating data export map." );

    // Build the data importer.
    d_data_export = Teuchos::rcp( new Tpetra::Export<GlobalOrdinal>(
				      d_export_map, d_import_map ) );
    testPostcondition( d_data_export != Teuchos::null,
		       "Error creating data importer." );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply the evaluation.
 */
template<class Mesh, class CoordinateField>
template<class SourceField, class TargetField>
void ConsistentEvaluation<Mesh,CoordinateField>::apply( 
    const Teuchos::RCP< FieldEvaluator<Mesh,SourceField> >& source_evaluator,
    TargetField& target_space )
{
    typedef FieldTraits<SourceField> SFT;
    typedef FieldTraits<TargetField> TFT;
    SourceField evaluated_field = 
	source_evaluator->evaluate( Teuchos::arcpFromArray( d_source_elements ), 
				    Teuchos::arcpFromArray( d_target_coords ) );
    testPrecondition( SFT::dim( evaluated_field ) == TFT::dim( target_space ),
		      "Source field dimension != target field dimension." );

    Teuchos::ArrayRCP<typename SFT::value_type> source_field_view;
    if ( SFT::size( evaluated_field ) == 0 )
    {
	source_field_view = 
	    Teuchos::ArrayRCP<typename SFT::value_type>( 0, 0.0 );
    }
    else
    {
	source_field_view = 
	    FieldTools<SourceField>::nonConstView( evaluated_field );
    }

    GlobalOrdinal source_size = SFT::size( evaluated_field ) /
				SFT::dim( evaluated_field );

    Teuchos::RCP< Tpetra::MultiVector<typename SFT::value_type,
				      GlobalOrdinal> > source_vector = 
	createMultiVectorFromView( d_export_map, 
				   source_field_view,
				   source_size,
				   SFT::dim( evaluated_field ) );

    Teuchos::ArrayRCP<typename TFT::value_type> target_field_view;
    if ( TFT::size( target_space ) == 0 )
    {
	target_field_view = 
	    Teuchos::ArrayRCP<typename TFT::value_type>( 0, 0.0 );
    }
    else
    {
	target_field_view =
	    FieldTools<TargetField>::nonConstView( target_space );
    }
    
    GlobalOrdinal target_size = TFT::size( target_space ) /
				TFT::dim( target_space );

    Teuchos::RCP< Tpetra::MultiVector<typename TFT::value_type,
				      GlobalOrdinal> > target_vector =	
	createMultiVectorFromView( d_import_map, 
				   target_field_view,
				   target_size,
				   TFT::dim( target_space ) );

    target_vector->doExport( *source_vector, *d_data_export, Tpetra::INSERT );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute globally unique ordinals for the points in the coordinate
 * field. 
 */
template<class Mesh, class CoordinateField>
Teuchos::Array<
    typename ConsistentEvaluation<Mesh,CoordinateField>::GlobalOrdinal>
ConsistentEvaluation<Mesh,CoordinateField>::computePointOrdinals(
    const CoordinateField& coordinate_field )
{
    int comm_rank = d_comm->getRank();
    int point_dim = CFT::dim( coordinate_field );
    GlobalOrdinal local_size = 
	std::distance( CFT::begin( coordinate_field ),
		       CFT::end( coordinate_field ) ) / point_dim;

    GlobalOrdinal global_size;
    Teuchos::reduceAll<int,GlobalOrdinal>( *d_comm,
					   Teuchos::REDUCE_MAX,
					   1,
					   &local_size,
					   &global_size );

    Teuchos::Array<GlobalOrdinal> point_ordinals( local_size );
    for ( int n = 0; n < local_size; ++n )
    {
	point_ordinals[n] = comm_rank*global_size + n;
    }
    return point_ordinals;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_CONSISTENTEVALUATION_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_ConsistentEvaluation_def.hpp
//---------------------------------------------------------------------------//

