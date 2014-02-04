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
 * \file DTK_IntegralAssemblyMap_def.hpp
 * \author Stuart R. Slattery
 * \brief Integral assembly map definition.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_INTEGRALASSEMBLYMAP_DEF_HPP
#define DTK_INTEGRALASSEMBLYMAP_DEF_HPP

#include <algorithm>
#include <limits>
#include <set>

#include "DTK_FieldTools.hpp"
#include "DTK_FieldTraits.hpp"
#include "DTK_DBC.hpp"
#include "DTK_Rendezvous.hpp"
#include "DTK_MeshTools.hpp"
#include "DTK_BoundingBox.hpp"

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_Import.hpp>
#include <Tpetra_Distributor.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>

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
 * \param geometric_tolerance Tolerance used for element vertex-in-geometry
 * checks.
 *
 * \param all_vertices_for_inclusion Flag for element-in-geometry
 * inclusion. If set to true, all of an element's vertices are required to
 * reside within a geometry within the geometric tolerance in order to be
 * considered a member of that geometry's conformal mesh. If set to false,
 * only one of an element's vertices must be contained within the geometric
 * tolerance of the geometry in order to be considered a member of that
 * geometry's conformal mesh.
 */
template<class Mesh, class Geometry>
IntegralAssemblyMap<Mesh,Geometry>::IntegralAssemblyMap(
    const RCP_Comm& comm, const int dimension, 
    const double geometric_tolerance, bool all_vertices_for_inclusion
 )
    : d_comm( comm )
    , d_dimension( dimension )
    , d_geometric_tolerance( geometric_tolerance )
    , d_all_vertices_for_inclusion( all_vertices_for_inclusion )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Mesh, class Geometry>
IntegralAssemblyMap<Mesh,Geometry>::~IntegralAssemblyMap()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Generate the integral map.
 *
 * \param source_mesh_manager Source mesh in the integral assembly problem. A
 * null RCP is a valid argument. This will be the case when a mesh manager is
 * only constructed on a subset of the processes that the integral assembly
 * map is constructed over. Note that the source mesh must exist only on
 * processes that reside within the IntegralAssemblyMap communicator.
 *
 * \param source_mesh_measure Measure interface object for getting mesh
 * element measures. A null RCP is a valid argument. This will be the case
 * when a mesh is only constructed on a subset of the processes that the
 * integral assembly map is constructed over.
 *
 * \param target_geometry_manager Target geometry in the integral assembly
 * problem. A null RCP is a valid argument. This will be the case when a
 * geometry manager is only constructed on a subset of the processes that the
 * integral assembly map is constructed over. Note that the target geometry
 * must exist only on processes that reside with the IntegralAssemblyMap
 * communicator. Geometry that exists outside this communication space will
 * not be considered in the mapping.
 */
template<class Mesh, class Geometry>
void IntegralAssemblyMap<Mesh,Geometry>::setup( 
    const RCP_MeshManager& source_mesh_manager, 
    const RCP_ElementMeasure& source_mesh_measure,
    const RCP_GeometryManager& target_geometry_manager )
{
    // Create existence values for the managers.
    bool source_exists = true;
    if ( source_mesh_manager.is_null() ) source_exists = false;
    bool target_exists = true;
    if ( target_geometry_manager.is_null() ) target_exists = false;
    d_comm->barrier();

    // Create local to global process indexers for the managers.
    RCP_Comm source_comm;
    if ( source_exists )
    {
	source_comm = source_mesh_manager->comm();
    }
    RCP_Comm target_comm;
    if ( target_exists )
    {
	target_comm = target_geometry_manager->comm();
    }
    d_comm->barrier();
    d_source_indexer = CommIndexer( d_comm, source_comm );
    d_target_indexer = CommIndexer( d_comm, target_comm );

    // Check the source and target dimensions for consistency.
    if ( source_exists )
    {
	DTK_REQUIRE( source_mesh_manager->dim() == d_dimension );
    }
    d_comm->barrier();

    if ( target_exists )
    {
	DTK_REQUIRE( target_geometry_manager->dim() == d_dimension );
    }
    d_comm->barrier();

    // Compute a unique global ordinal for each geometric object.
    Teuchos::Array<GlobalOrdinal> geometry_ordinals;
    computeGeometryOrdinals( target_geometry_manager, geometry_ordinals );

    // Build a rendezvous decomposition with the source mesh.
    BoundingBox global_box( -Teuchos::ScalarTraits<double>::rmax(),
			    -Teuchos::ScalarTraits<double>::rmax(),
			    -Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax(),
			    Teuchos::ScalarTraits<double>::rmax() );
    Rendezvous<Mesh> rendezvous( d_comm, d_dimension, global_box );
    rendezvous.build( source_mesh_manager );

    // Get the target geometries and their bounding boxes.
    Teuchos::ArrayRCP<Geometry> target_geometry(0);
    Teuchos::Array<BoundingBox> target_boxes(0);
    if ( target_exists )
    {
	target_boxes = target_geometry_manager->boundingBoxes();
	target_geometry = target_geometry_manager->geometry();
    }
    d_comm->barrier();

    // Allocate space for geometry data.
    d_integral_elements.resize( target_geometry.size() );
    d_geometry_measures.resize( target_geometry.size() );

    // Determine the rendezvous destination procs for the target geometries.
    Teuchos::Array<Teuchos::Array<int> > box_procs = 
	rendezvous.procsContainingBoxes( target_boxes );
    target_boxes.clear();

    // Unroll the rendezvous procs, target geometries, and ordinals to
    // populate the distrubutor.
    Teuchos::Array<int> rendezvous_procs;
    Teuchos::Array<Teuchos::Array<int> >::const_iterator box_iterator;
    Teuchos::Array<int>::const_iterator proc_iterator;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
	geometry_ordinal_iterator;
    Teuchos::Array<GlobalOrdinal> target_ordinals;
    Teuchos::Array<Geometry> target_geom_to_send;
    typename Teuchos::ArrayRCP<Geometry>::const_iterator 
	target_geometry_iterator;
    for ( box_iterator = box_procs.begin(), 
target_geometry_iterator = target_geometry.begin(),
geometry_ordinal_iterator = geometry_ordinals.begin();
	  box_iterator != box_procs.end();
	  ++box_iterator, ++target_geometry_iterator, 
			 ++geometry_ordinal_iterator )
    {
	for ( proc_iterator = box_iterator->begin();
	      proc_iterator != box_iterator->end();
	      ++proc_iterator )
	{
	    rendezvous_procs.push_back( *proc_iterator );
	    target_geom_to_send.push_back( *target_geometry_iterator );
	    target_ordinals.push_back( *geometry_ordinal_iterator );
	}
    }
    box_procs.clear();
    geometry_ordinals.clear();

    // Via an inverse communication, move the target geometries and their
    // ordinals to the rendezvous decomposition.
    Tpetra::Distributor target_to_rendezvous_distributor( d_comm );
    GlobalOrdinal num_rendezvous_geom = 
	target_to_rendezvous_distributor.createFromSends( rendezvous_procs() );
    rendezvous_procs.clear();

    Teuchos::ArrayView<const Geometry> target_geom_to_send_view =
	target_geom_to_send();
    Teuchos::Array<Geometry> 
	rendezvous_geometry( num_rendezvous_geom );
    target_to_rendezvous_distributor.doPostsAndWaits(
	target_geom_to_send_view, 1, rendezvous_geometry() );

    Teuchos::ArrayView<const GlobalOrdinal> target_ordinals_view =
	target_ordinals();
    Teuchos::Array<GlobalOrdinal> 
	rendezvous_geometry_ordinals( num_rendezvous_geom );
    target_to_rendezvous_distributor.doPostsAndWaits(
	target_ordinals_view, 1, rendezvous_geometry_ordinals() );

    // Extract the geometry target procs.    
    Teuchos::ArrayView<const int> from_images = 
	target_to_rendezvous_distributor.getImagesFrom();
    Teuchos::ArrayView<const std::size_t> from_lengths = 
	target_to_rendezvous_distributor.getLengthsFrom();
    Teuchos::Array<int> geom_target_procs;
    for ( int i = 0; i < Teuchos::as<int>(from_images.size()); ++i )
    {
	for ( std::size_t j = 0; j < from_lengths[i]; ++j )
	{
	    geom_target_procs.push_back( from_images[i] );
	}
    }
    DTK_CHECK( Teuchos::as<GlobalOrdinal>(geom_target_procs.size())
		   == num_rendezvous_geom );

    // Get the rendezvous source mesh elements that are in the rendezvous
    // target geometry. This is really expensive and we should rather think of
    // a way to logarithmically use the geometry bounding boxes for searching.
    Teuchos::Array<Teuchos::Array<GlobalOrdinal> > in_geom_elements;
    rendezvous.elementsInGeometry( rendezvous_geometry, in_geom_elements,
				   d_geometric_tolerance, 
				   d_all_vertices_for_inclusion );

    // Unroll the rendezvous elements and add the global target ordinal they
    // exist within.
    Teuchos::Array<GlobalOrdinal> rendezvous_elements;
    Teuchos::Array<GlobalOrdinal> in_geom_ordinals;
    Teuchos::Array<int> in_geom_target_procs;
    typename Teuchos::Array<Teuchos::Array<GlobalOrdinal> >::const_iterator 
	in_geom_elements_iterator;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator geom_element_iterator;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
	rendezvous_geom_ordinal_iterator;
    Teuchos::Array<int>::const_iterator geom_target_procs_iterator;
    for ( in_geom_elements_iterator = in_geom_elements.begin(),
   rendezvous_geom_ordinal_iterator = rendezvous_geometry_ordinals.begin(),
	 geom_target_procs_iterator = geom_target_procs.begin();
	  in_geom_elements_iterator != in_geom_elements.end();
	  ++in_geom_elements_iterator, ++rendezvous_geom_ordinal_iterator,
				      ++geom_target_procs_iterator )
    {
	for ( geom_element_iterator = in_geom_elements_iterator->begin();
	      geom_element_iterator != in_geom_elements_iterator->end();
	      ++geom_element_iterator )
	{
	    rendezvous_elements.push_back( *geom_element_iterator );
	    in_geom_ordinals.push_back( *rendezvous_geom_ordinal_iterator );
	    in_geom_target_procs.push_back( *geom_target_procs_iterator );
	}
    }
    in_geom_elements.clear();
    rendezvous_geometry_ordinals.clear();
    geom_target_procs.clear();

    // Communicate back to the target the elements that will construct the
    // integral for each geometry.
    Tpetra::Distributor rendezvous_to_target_distributor( d_comm );
    GlobalOrdinal num_target_elements = 
	rendezvous_to_target_distributor.createFromSends( 
	    in_geom_target_procs() );
    in_geom_target_procs.clear();

    Teuchos::ArrayView<const GlobalOrdinal> rendezvous_elements_view =
	rendezvous_elements();
    Teuchos::Array<GlobalOrdinal> mapped_target_elements( num_target_elements );
    rendezvous_to_target_distributor.doPostsAndWaits(
	rendezvous_elements_view, 1, mapped_target_elements() );

    Teuchos::ArrayView<const GlobalOrdinal> in_geom_ordinals_view =
	in_geom_ordinals();
    Teuchos::Array<GlobalOrdinal> mapped_target_ordinals( num_target_elements );
    rendezvous_to_target_distributor.doPostsAndWaits(
	in_geom_ordinals_view, 1, mapped_target_ordinals() );

    // Generate the list of integral elements.
    Teuchos::Array<std::set<GlobalOrdinal> > 
	integral_elements( target_geometry.size() );
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
	mapped_target_elements_it;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
	mapped_target_ordinals_it;
    for ( mapped_target_elements_it = mapped_target_elements.begin(),
	  mapped_target_ordinals_it = mapped_target_ordinals.begin();
	  mapped_target_elements_it != mapped_target_elements.end();
	  ++mapped_target_elements_it, ++mapped_target_ordinals_it )
    {
	DTK_CHECK( d_target_g2l.find(*mapped_target_ordinals_it) != 
		       d_target_g2l.end() );

	integral_elements[ 
	    d_target_g2l.find(*mapped_target_ordinals_it)->second ].
	    insert( *mapped_target_elements_it );
    }
    mapped_target_ordinals.clear();

    // Build a unique list of element ordinals in the target decomposition.
    typename Teuchos::Array<GlobalOrdinal>::iterator mapped_element_bound;
    std::sort( mapped_target_elements.begin(), mapped_target_elements.end() );
    mapped_element_bound = std::unique( mapped_target_elements.begin(),
					mapped_target_elements.end() );
    mapped_target_elements.resize( 
	std::distance( mapped_target_elements.begin(),
		       mapped_element_bound ) );

    // Build the element local to global map.
    std::tr1::unordered_map<GlobalOrdinal,GlobalOrdinal> element_g2l;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
	mapped_target_elements_begin = mapped_target_elements.begin();
    for ( mapped_target_elements_it = mapped_target_elements.begin();
	  mapped_target_elements_it != mapped_target_elements.end();
	  ++mapped_target_elements_it )
    {
	element_g2l[ *mapped_target_elements_it ] =
	    std::distance( mapped_target_elements_begin,
			   mapped_target_elements_it );
    }

    // Replace the global element ordinals in the integral element sets with
    // local ordinals.
    typename Teuchos::Array<std::set<GlobalOrdinal> >::iterator 
	integral_set_iterator;
    typename Teuchos::Array<Teuchos::Array<GlobalOrdinal> >::iterator 
	integral_elements_iterator;
    typename std::set<GlobalOrdinal>::iterator set_iterator;
    for ( integral_set_iterator = integral_elements.begin(),
     integral_elements_iterator = d_integral_elements.begin();
	  integral_set_iterator != integral_elements.end();
	  ++integral_set_iterator, ++integral_elements_iterator )
    {
	for ( set_iterator = integral_set_iterator->begin();
	      set_iterator != integral_set_iterator->end();
	      ++set_iterator )
	{
	    DTK_CHECK( element_g2l.find( *set_iterator ) != 
			   element_g2l.end() );

	    integral_elements_iterator->push_back( 
		element_g2l.find( *set_iterator )->second );
	}
    }
    integral_elements.clear();

    // Build the target map from the unique element ordinals.
    Teuchos::ArrayView<const GlobalOrdinal> mapped_target_elements_view = 
	mapped_target_elements();
    d_target_map = Tpetra::createNonContigMap<int,GlobalOrdinal>(
	mapped_target_elements_view, d_comm );
    DTK_ENSURE( !d_target_map.is_null() );

    // Allocate space for the element measures and get rid of the target
    // elements.
    Teuchos::Array<double> 
	integral_element_measures( mapped_target_elements.size() );
    mapped_target_elements.clear();

    // Put the rendezvous elements into a unique list and get their source
    // procs to populate the source distributor.
    typename Teuchos::Array<GlobalOrdinal>::iterator rendezvous_element_bound;
    std::sort( rendezvous_elements.begin(), rendezvous_elements.end() );
    rendezvous_element_bound = std::unique( rendezvous_elements.begin(),
					    rendezvous_elements.end() );
    rendezvous_elements.resize( std::distance( rendezvous_elements.begin(),
					       rendezvous_element_bound ) );
    Teuchos::Array<int> rendezvous_element_source_procs =
	rendezvous.elementSourceProcs( rendezvous_elements );

    // Communicate back to the source the elements we need integrals and
    // measures for.
    Teuchos::ArrayView<const int> rendezvous_element_source_procs_view =
	rendezvous_element_source_procs();
    Tpetra::Distributor rendezvous_to_source_distributor( d_comm );
    GlobalOrdinal num_source_elements = 
	rendezvous_to_source_distributor.createFromSends(
	    rendezvous_element_source_procs_view );

    d_source_elements.resize( num_source_elements );
    rendezvous_elements_view = rendezvous_elements();
    rendezvous_to_source_distributor.doPostsAndWaits( 
	rendezvous_elements_view, 1, d_source_elements() );

    // Build a unique list of the elements.
    typename Teuchos::Array<GlobalOrdinal>::iterator source_element_bound;
    std::sort( d_source_elements.begin(), d_source_elements.end() );
    source_element_bound = std::unique( d_source_elements.begin(),
					d_source_elements.end() );
    d_source_elements.resize( std::distance( d_source_elements.begin(),
					     source_element_bound ) );

    // Build the source map from the element ordinals.
    Teuchos::ArrayView<const GlobalOrdinal> source_elements_view =
	d_source_elements();
    d_source_map = Tpetra::createNonContigMap<int,GlobalOrdinal>(
	source_elements_view, d_comm );
    DTK_ENSURE( !d_source_map.is_null() );

    // Build the source-to-target importer.
    d_source_to_target_importer = 
      Teuchos::rcp( new Tpetra::Import<int,GlobalOrdinal>(
			  d_source_map, d_target_map ) );
    DTK_ENSURE( !d_source_to_target_importer.is_null() );

    // Communicate the element measures from the source to the target.
    Teuchos::Array<double> source_measures(0);
    if ( source_exists )
    {
	source_measures = source_mesh_measure->measure( 
	    Teuchos::arcpFromArray( d_source_elements ) );
    }
    Teuchos::RCP<Tpetra::Vector<double,int,GlobalOrdinal> > source_vector = 
	Tpetra::createVectorFromView( 
	    d_source_map, Teuchos::arcpFromArray( source_measures ) );
    Teuchos::RCP<Tpetra::Vector<double,int,GlobalOrdinal> > target_vector =
	Tpetra::createVectorFromView( 
	    d_target_map, 
	    Teuchos::arcpFromArray( integral_element_measures ) );
    target_vector->doImport( *source_vector, *d_source_to_target_importer,
			     Tpetra::INSERT );

    // Compute the local geometry measures from element measure sums. This
    // should approximate the true geometry measure.
    std::fill( d_geometry_measures.begin(), d_geometry_measures.end(), 0.0 );
    Teuchos::Array<double>::iterator geom_measure_iterator;
    typename Teuchos::Array<Teuchos::Array<GlobalOrdinal> >::const_iterator
	integral_iterator;
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
	integral_element_iterator;
    for ( geom_measure_iterator = d_geometry_measures.begin(),
	      integral_iterator = d_integral_elements.begin();
	  geom_measure_iterator != d_geometry_measures.end();
	  ++geom_measure_iterator, ++integral_iterator )
    {
	for ( integral_element_iterator = integral_iterator->begin();
	      integral_element_iterator != integral_iterator->end();
	      ++integral_element_iterator )
	{
	    *geom_measure_iterator +=
		integral_element_measures[*integral_element_iterator];
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply the integral assembly map for a valid source field integrator
 * and target data space to the target geometry.
 *
 * \param source_integrator Function integrator used to apply the
 * mapping. This FieldIntegrator must be valid for the source mesh used to
 * generate the map.
 *
 * \param target_space_manager Target space into which the function
 * integrations will be written. Enough space must be allocated to hold
 * integrations in all geometries in all dimensions of the field.
 */
template<class Mesh, class Geometry>
template<class SourceField, class TargetField>
void IntegralAssemblyMap<Mesh,Geometry>::apply( 
    const Teuchos::RCP<FieldIntegrator<Mesh,SourceField> >& source_integrator,
    Teuchos::RCP<FieldManager<TargetField> >& target_space_manager )
{
    typedef FieldTraits<SourceField> SFT;
    typedef FieldTraits<TargetField> TFT;

    // Set existence values for the source and target.
    bool source_exists = true;
    if ( source_integrator.is_null() ) source_exists = false;
    bool target_exists = true;
    if ( target_space_manager.is_null() ) target_exists = false;
    d_comm->barrier();

    // Integrate the source function at the target points and construct a view
    // of the function integrals.
    int source_dim;
    Teuchos::ArrayRCP<typename SFT::value_type> source_field_copy(0,0);
    if ( source_exists )
    {
	SourceField function_integrations = 
	    source_integrator->integrate( 
		Teuchos::arcpFromArray( d_source_elements ) );

	source_dim = SFT::dim( function_integrations );

	source_field_copy =    
	    FieldTools<SourceField>::copy( function_integrations );
    }
    d_comm->barrier();
    Teuchos::broadcast<int,int>( *d_comm, d_source_indexer.l2g(0),
				 Teuchos::Ptr<int>(&source_dim) );

    // Build a multivector for the function integrations.
    GlobalOrdinal source_size = source_field_copy.size() / source_dim;
    Teuchos::RCP<Tpetra::MultiVector<typename SFT::value_type,int,GlobalOrdinal> > 
	source_vector = Tpetra::createMultiVectorFromView( 
	    d_source_map, source_field_copy, source_size, source_dim );

    // Construct a view of the target space and fill it with zeros so that we
    // start the integral summations at zero.
    GlobalOrdinal integral_size = d_geometry_measures.size();
    int target_dim;
    DTK_REMEMBER( GlobalOrdinal target_size = 0 );
    Teuchos::ArrayRCP<typename TFT::value_type> target_field_view(0,0);
    if ( target_exists )
    {
	target_field_view = FieldTools<TargetField>::nonConstView( 
	    *target_space_manager->field() );

	target_dim = TFT::dim( *target_space_manager->field() );

	DTK_REMEMBER( target_size = target_dim * integral_size );

	FieldTools<TargetField>::putScalar( 
	    *target_space_manager->field(), 0.0 );
    }
    d_comm->barrier();
    Teuchos::broadcast<int,int>( *d_comm, d_target_indexer.l2g(0),
				 Teuchos::Ptr<int>(&target_dim) );

    // Check that the source and target have the same field dimension.
    DTK_REQUIRE( source_dim == target_dim );

    // Verify that the target space has the proper amount of memory allocated.
    DTK_REQUIRE( target_size == 
		      Teuchos::as<GlobalOrdinal>(target_field_view.size()) );

    // Build a multivector for the function integrations in the target
    // decomposition.
    Tpetra::MultiVector<typename TFT::value_type,int,GlobalOrdinal> 
	target_vector( d_target_map, target_dim );

    // Import the function integrations.
    target_vector.doImport( *source_vector, *d_source_to_target_importer,
			     Tpetra::INSERT );

    // Collapse the function integrations over the geometry and apply them to
    // the target space.
    typename Teuchos::Array<GlobalOrdinal>::const_iterator 
	integral_element_iterator;
    typename TFT::size_type target_index;
    Teuchos::ArrayRCP<const typename TFT::value_type> dim_element_integrals;
    for ( typename TFT::size_type n = 0; 
	  n < Teuchos::as<typename TFT::size_type>(integral_size); ++n )
    {
	for ( int d = 0; d < target_dim; ++d )
	{
	    dim_element_integrals = target_vector.getData( d );
	    target_index = n + integral_size*d;

	    for ( integral_element_iterator = d_integral_elements[n].begin();
		  integral_element_iterator != d_integral_elements[n].end();
		  ++integral_element_iterator )
	    {
		target_field_view[ target_index ] +=
		    dim_element_integrals[*integral_element_iterator];
	    }
	}
    }

    // Scale the results by the inverse geometry measure sums.
    for ( typename TFT::size_type n = 0; 
	  n < Teuchos::as<typename TFT::size_type>(integral_size); ++n )
    {
	for ( int d = 0; d < target_dim; ++d )
	{
	    target_index = n + integral_size*d;
	    target_field_view[ target_index ] /=
		Teuchos::as<typename TFT::value_type>(d_geometry_measures[n]);
	}
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute globally unique ordinals for the target points. Here an
 * invalid ordinal will be designated as the maximum value as specified by the
 * limits header for the ordinal type. We do this so that 0 may be a valid
 * ordinal.
 *
 * \param target_geometry_manager The geometries to compute global ordinals for.
 *
 * \param target_ordinals The computed globally unique ordinals for the target
 * geometries. 
 */
template<class Mesh, class Geometry>
void IntegralAssemblyMap<Mesh,Geometry>::computeGeometryOrdinals(
    const RCP_GeometryManager& target_geometry_manager,
    Teuchos::Array<GlobalOrdinal>& target_ordinals )
{
    // Set an existence value for the target geometry.
    bool target_exists = true;
    if ( target_geometry_manager.is_null() ) target_exists = false;
    int comm_rank = d_comm->getRank();
    GlobalOrdinal local_size = 0;

    if ( target_exists )
    {
	local_size = target_geometry_manager->localNumGeometry();
    }
    d_comm->barrier();

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
	d_target_g2l[ target_ordinals[n] ] = n;
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

#endif // end DTK_INTEGRALASSEMBLYMAP_DEF_HPP

//---------------------------------------------------------------------------//
// end DTK_IntegralAssemblyMap_def.hpp
//---------------------------------------------------------------------------//
