//---------------------------------------------------------------------------//
/*!
 * \file DTK_ConsistentInterpolation.hpp
 * \author Stuart R. Slattery
 * \brief Consistent interpolation mapping declaration.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_CONSISTENTINTERPOLATION_HPP
#define DTK_CONSISTENTINTERPOLATION_HPP

#include <set>

#include <DTK_Exception.hpp>
#include <DTK_Rendezvous.hpp>
#include <DTK_MeshTools.hpp>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ENull.hpp>

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
    BoundingBox coord_box = FieldTools<CoordinateField>::coordGlobalBoundingBox(
	coordinate_field, d_comm );

    // Intersect the boxes to get the rendezvous bounding box.
    BoundingBox rendezvous_box = buildRendezvousBox( 
	mesh, mesh_box, coordinate_field, coord_box );

    // Build a rendezvous decomposition.
    Rendezvous<Mesh> rendezvous( d_comm, rendezvous_box );

    // Compute a unique global ordinal for each point in the coordinate field.
    Teuchos::Array<global_ordinal_type> point_ordinals = 
	computePointOrdinals( coordinate_field );

    // Build the import map from the global ordinals.
    Teuchos::ArrayView<const global_ordinal_type> import_ordinal_view =
	point_ordinals();
    d_import_map = Tpetra::createNonContigMap<global_ordinal_type>(
	import_ordinal_view, d_comm );
    testPostcondition( d_import_map != Teuchos::null,
		       "Error creating data import map." );
    
    // Determine the rendezvous destination of each point in the coordinate
    // field.
    std::size_t coord_dim = CFT::dim( coordinate_field );
    CFT::size_type num_coords = CFT::size( coordinate_field );
    Teuchos::ArrayRCP coords_view( &*CFT::begin( coordinate_field ), 
				   num_coords );
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

    // Setup Target to rendezvous communication.
    Teuchos::ArrayView<const global_ordinal_type> rendezvous_points_view =
	rendezvous_points();
    RCP_TpetraMap rendezvous_import_map = 
	Tpetra::createNonContigMap<global_ordinal_type>(
	    rendezvous_points_view, d_comm );
    Tpetra::Export point_exporter( d_import_map, rendezvous_import_map );
    
    // Move the target coordinates to the rendezvous decomposition.
    global_ordinal_type num_points = num_coords / coord_dim;
    Teuchos::RCP< Tpetra::MultiVector<double,global_ordinal_type> > 
	target_coords =	createMultiVectorFromView( d_import_map, coords_view,
						   num_points, coord_dim );
    Tpetra::MultiVector<double,global_ordinal_type> rendezvous_coords(
	rendezvous_import_map, coord_dim );
    rendezvous_coords.doExport( *target_coords, point_exporter, Tpetra::INSERT );

    // Search the rendezvous decomposition with the points.
    Teuchos::Array<global_ordinal_type> rendezvous_elements =
	rendezvous.getElements( rendezvous_coords.get1dView() );

    // Setup rendezvous to source communication.
    global_ordinal_type num_elements = std::distance( MT::elementsBegin( mesh ), 
						      MT::elementsEnd( mesh ) );
    Teuchos::ArrayView<const global_ordinal_type> mesh_element_view(
	&*MT::elementsBegin( mesh ), num_elements );
    RCP_TpetraMap mesh_element_map = 
	Tpetra::createNonContigMap<global_ordinal_type>( 
	    mesh_element_view, d_comm );

    Teuchos::ArrayView<const global_ordinal_type> rendezvous_elements_view =
	rendezvous_elements();
    RCP_TpetraMap rendezvous_element_map = 
	Tpetra::createNonContigMap<global_ordinal_type>( 
	    rendezvous_elements_view, d_comm );

    Tpetra::Import<global_ordinal_type> source_importer( rendezvous_element_map,
							 mesh_element_map );

    // Send the point coordinates to the source decomposition.

    // Send the point ordinals to the source decomposition.
    
    // Build the data export map from the coordinate ordinals.

    // Build the data importer.
    d_data_importer = Teuchos::rcp( new Tpetra::Import<global_ordinal_type>(
					d_export_map, d_import_map ) );
    testPostcondition( d_data_importer != Teuchos::null,
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
	source_evaluator->evaluate( d_source_elements, d_target_coords );

    testPrecondition( SFT::dim( evaluated_field ) == TFT::dim( target_space ),
		      "Source field dimension != target field dimension." );

    Teuchos::ArrayRCP<SFT::value_type> 
	source_field_view( &*SFT::begin( evaluated_field ),
			   SFT::size( evaluated_field ) );

    global_ordinal_type source_size = SFT::size( evaluated_field ) /
				      SFT::dim( evaluated_field );

    Teuchos::RCP< Tpetra::MultiVector<SFT::value_type> > source_vector =
	createMultiVectorFromView( d_export_map, 
				   source_field_view,
				   source_size,
				   SFT::dim( evaluated_field ) );

    Teuchos::ArrayRCP<TFT::value_type> target_field_view(
	&*TFT::begin( target_space ),
	TFT::size( target_space ) );
    
    global_ordinal_type target_size = TFT::size( target_space ) /
				      FFT::dim( target_space );

    Teuchos::RCP< Tpetra::MultiVector<TFT::value_type> > target_vector =
	createMultiVectorFromView( d_import_map, 
				   target_field_view,
				   target_size,
				   TFT::dim( target_space ) );

    target_vector->doImport( *source_vector, *d_data_importer, Tpetra::INSERT );
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
    CFT::size_type num_points = 
	std::distance( CFT::begin( coordinate_field ),
		       CFT::end( coordinate_field ) ) / point_dim;
    CFT::const_iterator point_coords = CFT::begin( coordinate_field );
    double point[3];
    for ( CFT::size_type n = 0; n < num_points; ++n )
    {
	for ( std::size_t d = 0; d < point_dim; ++d )
	{
	    point[d] = point_coords[d*num_points + n];
	}
	for ( std::size_t d = point_dim; d < 3; ++d )
	{
	    point[d] = 0.0;
	}

	if ( pointInBox( point ) )
	{
	    x_values.insert( point[0] );
	    y_values.insert( point[1] );
	    z_values.insert( point[2] );
	}
    }

    // Get the mesh nodes in the coord box.
    std::size_t node_dim = MT::nodeDim( mesh );
    MT::size_type num_nodes = 
	std::distance( MT::begin( mesh ), MT::end( mesh ) ) / node_dim;
    MT::const_iterator node_coords = MT::coordsBegin( mesh );
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

	if ( nodeInBox( node ) )
	{
	    x_values.insert( node[0] );
	    y_values.insert( node[1] );
	    z_values.insert( node[2] );
	}
    }

    // Compute the global bounding box for the collected coordinates.
    double global_x_min, global_y_min, global_z_min;
    double global_x_max, global_y_max, global_z_max;

    Teuchos::reduceAll<int,double>( *d_comm, 
				    Teuchos::REDUCE_MIN,
				    x_values.size(),
				    &x_values[0],
				    &global_x_min );

    Teuchos::reduceAll<int,double>( *d_comm, 
				    Teuchos::REDUCE_MIN,
				    y_values.size(),
				    &y_values[0],
				    &global_y_min );

    Teuchos::reduceAll<int,double>( *d_comm, 
				    Teuchos::REDUCE_MIN,
				    z_values.size(),
				    &z_values[0],
				    &global_z_min );

    Teuchos::reduceAll<int,double>( *d_comm, 
				    Teuchos::REDUCE_MAX,
				    x_values.size(),
				    &x_values[0],
				    &global_x_max );

    Teuchos::reduceAll<int,double>( *d_comm, 
				    Teuchos::REDUCE_MAX,
				    y_values.size(),
				    &y_values[0],
				    &global_y_max );

    Teuchos::reduceAll<int,double>( *d_comm, 
				    Teuchos::REDUCE_MAX,
				    z_values.size(),
				    &z_values[0],
				    &global_z_max );

    return BoundingBox( global_x_min, global_y_min, global_z_min,
			global_x_max, global_y_max, global_z_max );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute globally unique ordinals for the points in the coordinate
 * field. 
 */
template<class Mesh, class CoordinateField>
Teuchos::Array<global_ordinal_type> 
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
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_CONSISTENTINTERPOLATION_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_ConsistentInterpolation_def.hpp
//---------------------------------------------------------------------------//

