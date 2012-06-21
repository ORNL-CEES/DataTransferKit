//---------------------------------------------------------------------------//
/*!
 * \file DTK_ConsistentInterpolation.hpp
 * \author Stuart R. Slattery
 * \brief Consistent interpolation mapping declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CONSISTENTINTERPOLATION_DEF_HPP
#define DTK_CONSISTENTINTERPOLATION_DEF_HPP

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
ConsistentInterpolation<Mesh,CoordinateField>::ConsistentInterpolation( 
    const RCP_Comm& comm )
    : d_comm( comm )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Mesh, class CoordinateField>
ConsistentInterpolation<Mesh,CoordinateField>::~ConsistentInterpolation()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Setup for interpolation.
 */
template<class Mesh, class CoordinateField>
void ConsistentInterpolation<Mesh,CoordinateField>::setup( 
    const Mesh& mesh, const CoordinateField& coordinate_field )
{
    // Get the global bounding box for the mesh.
    BoundingBox mesh_box = MeshTools<Mesh>::globalBoundingBox( mesh, d_comm );

    // Get the global bounding box for the coordinate field.
    BoundingBox coord_box = 
	FieldTools<CoordinateField>::coordGlobalBoundingBox(
	    coordinate_field, d_comm );

    // Intersect the boxes to get the rendezvous bounding box.
    BoundingBox rendezvous_box = buildRendezvousBox( 
	mesh, mesh_box, coordinate_field, coord_box );

    // Build a rendezvous decomposition with the source mesh.
    Rendezvous<Mesh> rendezvous( d_comm, rendezvous_box );
    rendezvous.build( mesh );

    // Compute a unique global ordinal for each point in the coordinate field.
    Teuchos::Array<global_ordinal_type> point_ordinals = 
	computePointOrdinals( coordinate_field );

    // Build the data import map from the point global ordinals.
    Teuchos::ArrayView<const global_ordinal_type> import_ordinal_view =
	point_ordinals();
    d_import_map = Tpetra::createNonContigMap<global_ordinal_type>(
	import_ordinal_view, d_comm );
    testPostcondition( d_import_map != Teuchos::null,
		       "Error creating data import map." );

    // Determine the rendezvous destination proc of each point in the
    // coordinate field.
    std::size_t coord_dim = CFT::dim( coordinate_field );
    typename CFT::size_type num_coords = CFT::size( coordinate_field );
    Teuchos::ArrayRCP<double> coords_view;
    if ( num_coords == 0 )
    {
	coords_view = Teuchos::ArrayRCP<double>( 0, 0.0 );
    }
    else
    {
	 coords_view = Teuchos::ArrayRCP<double>( 
	     (double*) &*CFT::begin( coordinate_field ), 0, num_coords, false );
    }
    Teuchos::Array<int> rendezvous_procs = 
	rendezvous.getRendezvousProcs( coords_view );

    // Via an inverse communication operation, move the global point ordinals
    // to the rendezvous decomposition.
    Tpetra::Distributor point_distributor( d_comm );
    global_ordinal_type num_rendezvous_points = 
	point_distributor.createFromSends( rendezvous_procs() );
    Teuchos::Array<global_ordinal_type> 
	rendezvous_points( num_rendezvous_points );
    point_distributor.doPostsAndWaits( 
	import_ordinal_view, 1, rendezvous_points() );

    // Setup Target-to-rendezvous communication.
    Teuchos::ArrayView<const global_ordinal_type> rendezvous_points_view =
	rendezvous_points();
    RCP_TpetraMap rendezvous_import_map = 
	Tpetra::createNonContigMap<global_ordinal_type>(
	    rendezvous_points_view, d_comm );
    Tpetra::Export<global_ordinal_type> 
	point_exporter( d_import_map, rendezvous_import_map );

    // Move the target coordinates to the rendezvous decomposition.
    global_ordinal_type num_points = num_coords / coord_dim;
    Teuchos::RCP< Tpetra::MultiVector<double,global_ordinal_type> > 
	target_coords =	createMultiVectorFromView( d_import_map, coords_view,
						   num_points, coord_dim );
    Tpetra::MultiVector<double,global_ordinal_type> rendezvous_coords( 
	rendezvous_import_map, coord_dim );
    rendezvous_coords.doExport( *target_coords, point_exporter, 
				Tpetra::INSERT );

    // Search the rendezvous decomposition with the target points to get the
    // source elements that contain them.
    Teuchos::Array<global_ordinal_type> rendezvous_elements =
	rendezvous.getElements( rendezvous_coords.get1dViewNonConst() );

    // Build a unique list of rendezvous elements.
    Teuchos::Array<global_ordinal_type> rendezvous_element_set =
	rendezvous_elements;
    std::sort( rendezvous_element_set.begin(), rendezvous_element_set.end() );
    typename Teuchos::Array<global_ordinal_type>::iterator unique_bound =
	std::unique( rendezvous_element_set.begin(), 
		     rendezvous_element_set.end() );
    rendezvous_element_set.resize( unique_bound - rendezvous_element_set.begin() );

    // Setup source-to-rendezvous communication.
    Teuchos::ArrayView<const global_ordinal_type> rendezvous_elem_set_view =
	rendezvous_element_set();
    RCP_TpetraMap rendezvous_element_map = 
	Tpetra::createNonContigMap<global_ordinal_type>( 
	    rendezvous_elem_set_view, d_comm );

    global_ordinal_type num_elements = 
	std::distance( MT::elementsBegin( mesh ), MT::elementsEnd( mesh ) );
    Teuchos::ArrayView<const global_ordinal_type> mesh_element_view(
	&*MT::elementsBegin( mesh ), num_elements );
    RCP_TpetraMap mesh_element_map = 
	Tpetra::createNonContigMap<global_ordinal_type>( 
	    mesh_element_view, d_comm );

    Tpetra::Import<global_ordinal_type> source_importer( 
	mesh_element_map, rendezvous_element_map );

    // Elements send their source decomposition proc to the rendezvous
    // decomposition.
    Tpetra::MultiVector<int,global_ordinal_type> 
	source_procs( mesh_element_map, 1 );
    source_procs.putScalar( d_comm->getRank() );
    Tpetra::MultiVector<int,global_ordinal_type> 
	rendezvous_source_procs( rendezvous_element_map, 1 );
    rendezvous_source_procs.doImport( source_procs, source_importer, 
				      Tpetra::INSERT );

    // Set the destination procs for the rendezvous-to-source communication.
    Teuchos::ArrayRCP<const int> rendezvous_source_procs_view =
	rendezvous_source_procs.get1dView();
    global_ordinal_type num_rendezvous_elements = rendezvous_elements.size();
    global_ordinal_type dest_proc_idx;
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
    Teuchos::ArrayView<const global_ordinal_type> rendezvous_elements_view =
	rendezvous_elements();
    Tpetra::Distributor source_distributor( d_comm );
    global_ordinal_type num_source_elements = 
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
    Teuchos::Array<global_ordinal_type> source_points( num_source_elements );
    source_distributor.doPostsAndWaits( rendezvous_points_view, 1,
					source_points() );

    // Build the data export map from the coordinate ordinals as well as
    // populate the source element/target coordinate pairings.
    Teuchos::ArrayView<const global_ordinal_type> source_points_view =
	source_points();
    d_export_map = Tpetra::createNonContigMap<global_ordinal_type>(
	source_points_view, d_comm );
    testPostcondition( d_export_map != Teuchos::null,
		       "Error creating data export map." );

    // Build the data importer.
    d_data_export = Teuchos::rcp( new Tpetra::Export<global_ordinal_type>(
					d_export_map, d_import_map ) );
    testPostcondition( d_data_export != Teuchos::null,
		       "Error creating data importer." );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply the interpolation.
 */
template<class Mesh, class CoordinateField>
template<class SourceField, class TargetField>
void ConsistentInterpolation<Mesh,CoordinateField>::apply( 
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
	 source_field_view = Teuchos::ArrayRCP<typename SFT::value_type>(
	     &*SFT::begin( evaluated_field ), 0, 
	     SFT::size( evaluated_field ), false );
    }

    global_ordinal_type source_size = SFT::size( evaluated_field ) /
				      SFT::dim( evaluated_field );

    Teuchos::RCP< Tpetra::MultiVector<typename SFT::value_type,global_ordinal_type> > 
	source_vector = createMultiVectorFromView( d_export_map, 
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
	target_field_view = Teuchos::ArrayRCP<typename TFT::value_type>(
	    &*TFT::begin( target_space ), 0,
	    TFT::size( target_space ), false );
    }
    
    global_ordinal_type target_size = TFT::size( target_space ) /
				      TFT::dim( target_space );

    Teuchos::RCP< Tpetra::MultiVector<typename TFT::value_type,global_ordinal_type> > 
	target_vector =	createMultiVectorFromView( d_import_map, 
						   target_field_view,
						   target_size,
						   TFT::dim( target_space ) );

    target_vector->doExport( *source_vector, *d_data_export, Tpetra::INSERT );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the bounding box for the rendezvous decomposition.
 */
template<class Mesh, class CoordinateField>
BoundingBox ConsistentInterpolation<Mesh,CoordinateField>::buildRendezvousBox( 
    const Mesh& mesh, const BoundingBox& mesh_box,
    const CoordinateField& coordinate_field, const BoundingBox& coord_box )
{
    std::set<double> x_values, y_values, z_values;
    
    // Get the coords in the mesh box.
    std::size_t point_dim = CFT::dim( coordinate_field );
    typename CFT::size_type num_points = 
	std::distance( CFT::begin( coordinate_field ),
		       CFT::end( coordinate_field ) ) / point_dim;
    typename CFT::const_iterator point_coords = CFT::begin( coordinate_field );
    double point[3];
    for ( typename CFT::size_type n = 0; n < num_points; ++n )
    {
	for ( std::size_t d = 0; d < point_dim; ++d )
	{
	    point[d] = point_coords[d*num_points + n];
	}
	for ( std::size_t d = point_dim; d < 3; ++d )
	{
	    point[d] = 0.0;
	}

	if ( mesh_box.pointInBox( point ) )
	{
	    x_values.insert( point[0] );
	    y_values.insert( point[1] );
	    z_values.insert( point[2] );
	}
    }

    // Get the mesh nodes in the coord box.
    std::size_t node_dim = MT::nodeDim( mesh );
    typename MT::global_ordinal_type num_nodes = 
	std::distance( MT::nodesBegin( mesh ), MT::nodesEnd( mesh ) ) / node_dim;
    typename MT::const_coordinate_iterator node_coords = MT::coordsBegin( mesh );
    double node[3];
    for ( global_ordinal_type n = 0; n < num_nodes; ++n )
    {
	for ( std::size_t d = 0; d < node_dim; ++d )
	{
	    node[d] = node_coords[d*num_nodes + n];
	}
	for ( std::size_t d = node_dim; d < 3; ++d )
	{
	    node[d] = 0.0;
	}

	if ( coord_box.pointInBox( node ) )
	{
	    x_values.insert( node[0] );
	    y_values.insert( node[1] );
	    z_values.insert( node[2] );
	}
    }
    
    // Compute the local min and max values.
    double local_x_min = *x_values.begin();
    double local_y_min = *y_values.begin();
    double local_z_min = *z_values.begin();

    double local_x_max = *x_values.rbegin();
    double local_y_max = *y_values.rbegin();
    double local_z_max = *z_values.rbegin();

    x_values.clear();
    y_values.clear();
    z_values.clear();

    // Compute the global bounding box for the collected coordinates.
    double global_x_min = -Teuchos::ScalarTraits<double>::rmax();
    double global_y_min = -Teuchos::ScalarTraits<double>::rmax();
    double global_z_min = -Teuchos::ScalarTraits<double>::rmax();

    double global_x_max = Teuchos::ScalarTraits<double>::rmax();
    double global_y_max = Teuchos::ScalarTraits<double>::rmax();
    double global_z_max = Teuchos::ScalarTraits<double>::rmax();
    
    if ( point_dim > 0 )
    {
	Teuchos::reduceAll<int,double>( *d_comm, 
					Teuchos::REDUCE_MIN,
					1,
					&local_x_min,
					&global_x_min );

	Teuchos::reduceAll<int,double>( *d_comm, 
					Teuchos::REDUCE_MAX,
					1,
					&local_x_max,
					&global_x_max );
    }
    if ( point_dim > 1 )
    {
	Teuchos::reduceAll<int,double>( *d_comm, 
					Teuchos::REDUCE_MIN,
					1,
					&local_y_min,
					&global_y_min );

	Teuchos::reduceAll<int,double>( *d_comm, 
					Teuchos::REDUCE_MAX,
					1,
					&local_y_max,
					&global_y_max );
    }
    if ( point_dim > 2 )
    {
	Teuchos::reduceAll<int,double>( *d_comm, 
					Teuchos::REDUCE_MIN,
					1,
					&local_z_min,
					&global_z_min );

	Teuchos::reduceAll<int,double>( *d_comm, 
					Teuchos::REDUCE_MAX,
					1,
					&local_z_max,
					&global_z_max );
    }

    // return BoundingBox( global_x_min, global_y_min, global_z_min,
    // 			global_x_max, global_y_max, global_z_max );

    // Right now I'm returning the mesh box. Scaling study 1 made apparent a
    // bad error in the computation of the bounding box intersection computed
    // by the code above. I'll need to update the code to handle corner cases
    // like this. For now, this will handle most cases but at the cost of
    // repartitioning the entire mesh when only a subset will be needed.
    return mesh_box;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute globally unique ordinals for the points in the coordinate
 * field. 
 */
template<class Mesh, class CoordinateField>
Teuchos::Array<
typename ConsistentInterpolation<Mesh,CoordinateField>::global_ordinal_type>
ConsistentInterpolation<Mesh,CoordinateField>::computePointOrdinals(
    const CoordinateField& coordinate_field )
{
    int comm_rank = d_comm->getRank();
    int point_dim = CFT::dim( coordinate_field );
    global_ordinal_type local_size = 
	std::distance( CFT::begin( coordinate_field ),
		       CFT::end( coordinate_field ) ) / point_dim;

    global_ordinal_type global_size;
    Teuchos::reduceAll<int,global_ordinal_type>( *d_comm,
					    Teuchos::REDUCE_MAX,
					    1,
					    &local_size,
					    &global_size );

    Teuchos::Array<global_ordinal_type> point_ordinals( local_size );
    for ( int n = 0; n < local_size; ++n )
    {
	point_ordinals[n] = comm_rank*global_size + n;
    }
    return point_ordinals;
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_CONSISTENTINTERPOLATION_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_ConsistentInterpolation_def.hpp
//---------------------------------------------------------------------------//

