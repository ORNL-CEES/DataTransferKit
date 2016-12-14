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
 * \file   DTK_MovingLeastSquareReconstructionOperator_impl.hpp
 * \author Stuart R. Slattery
 * \brief  Moving least square interpolator.
 */
//---------------------------------------------------------------------------//

#ifndef DTK_MOVINGLEASTSQUARERECONSTRUCTIONOPERATOR_IMPL_HPP
#define DTK_MOVINGLEASTSQUARERECONSTRUCTIONOPERATOR_IMPL_HPP

#include "DTK_BasicEntityPredicates.hpp"
#include "DTK_CenterDistributor.hpp"
#include "DTK_DBC.hpp"
#include "DTK_LocalMLSProblem.hpp"
#include "DTK_MovingLeastSquareReconstructionOperator.hpp"
#include "DTK_PredicateComposition.hpp"
#include "DTK_SplineInterpolationPairing.hpp"

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Ptr.hpp>

#include <Tpetra_MultiVector.hpp>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Constructor.
template <class Basis, int DIM>
MovingLeastSquareReconstructionOperator<Basis, DIM>::
    MovingLeastSquareReconstructionOperator(
        const Teuchos::RCP<const TpetraMap> &domain_map,
        const Teuchos::RCP<const TpetraMap> &range_map,
        const Teuchos::ParameterList &parameters )
    : Base( domain_map, range_map )
    , d_use_knn( false )
    , d_knn( 0 )
    , d_radius( 0.0 )
    , d_domain_entity_dim( 0 )
    , d_range_entity_dim( 0 )
{
    // Determine if we are doing kNN search or radius search.
    if ( parameters.isParameter( "Type of Search" ) )
    {
        if ( "Radius" == parameters.get<std::string>( "Type of Search" ) )
        {
            d_use_knn = false;
        }
        else if ( "Nearest Neighbor" ==
                  parameters.get<std::string>( "Type of Search" ) )
        {
            d_use_knn = true;
        }
        else
        {
            // Otherwise we got an invalid search type.
            DTK_INSIST( false );
        }
    }

    // If we are doing kNN support get the number of neighbors.
    if ( d_use_knn )
    {
        DTK_REQUIRE( parameters.isParameter( "Num Neighbors" ) );
        d_knn = parameters.get<int>( "Num Neighbors" );
    }

    // Otherwise we are doing the radius search so get the basis radius.
    else
    {
        DTK_REQUIRE( parameters.isParameter( "RBF Radius" ) );
        d_radius = parameters.get<double>( "RBF Radius" );
    }

    // Get the topological dimension of the domain and range entities. This
    // map will use their centroids for the point cloud.
    if ( parameters.isParameter( "Domain Entity Dimension" ) )
    {
        d_domain_entity_dim = parameters.get<int>( "Domain Entity Dimension" );
    }
    if ( parameters.isParameter( "Range Entity Dimension" ) )
    {
        d_range_entity_dim = parameters.get<int>( "Range Entity Dimension" );
    }
}

//---------------------------------------------------------------------------//
// Setup the map operator.
template <class Basis, int DIM>
void MovingLeastSquareReconstructionOperator<Basis, DIM>::setupImpl(
    const Teuchos::RCP<FunctionSpace> &domain_space,
    const Teuchos::RCP<FunctionSpace> &range_space )
{
    DTK_REQUIRE( Teuchos::nonnull( domain_space ) );
    DTK_REQUIRE( Teuchos::nonnull( range_space ) );

    // Extract the Support maps.
    const Teuchos::RCP<const typename Base::TpetraMap> domain_map =
        this->getDomainMap();
    const Teuchos::RCP<const typename Base::TpetraMap> range_map =
        this->getRangeMap();

    // Get the parallel communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> comm = domain_map->getComm();

    // Extract the source nodes and their ids.
    Teuchos::ArrayRCP<double> source_centers;
    Teuchos::ArrayRCP<GO> source_support_ids;
    getNodeCoordsAndIds( domain_space, d_domain_entity_dim, source_centers,
                         source_support_ids );

    // Extract the target nodes and their ids.
    Teuchos::ArrayRCP<double> target_centers;
    Teuchos::ArrayRCP<GO> target_support_ids;
    getNodeCoordsAndIds( range_space, d_range_entity_dim, target_centers,
                         target_support_ids );

    // Calculate an approximate neighborhood distance for the local target
    // centers. If using kNN, compute an approximation. If doing a radial
    // search, use the radius. We will use these distances to expand the local
    // bounding box to ensure we find all of our neighbors in parallel.
    double target_proximity = 0.0;
    if ( d_use_knn )
    {
        // Get the local bounding box.
        Teuchos::Tuple<double, 6> local_box;
        range_space->entitySet()->localBoundingBox( local_box );

        // Calculate the largest span of the cardinal directions.
        target_proximity = local_box[3] - local_box[0];
        for ( int d = 1; d < DIM; ++d )
        {
            target_proximity =
                std::max( target_proximity, local_box[d + 3] - local_box[d] );
        }

        // Take the proximity to be 10% of the largest distance.
        target_proximity *= 0.1;
    }
    else
    {
        target_proximity = d_radius;
    }

    // Gather the source centers that are in the proximity of the target
    // centers on this proc.
    Teuchos::Array<double> dist_sources;
    CenterDistributor<DIM> distributor( comm, source_centers(),
                                        target_centers(), target_proximity,
                                        dist_sources );

    // Gather the global ids of the source centers that are within the proximity
    // of
    // the target centers on this proc.
    Teuchos::Array<GO> dist_source_support_ids( distributor.getNumImports() );
    Teuchos::ArrayView<const GO> source_support_ids_view = source_support_ids();
    distributor.distribute( source_support_ids_view,
                            dist_source_support_ids() );

    // Build the source/target pairings.
    SplineInterpolationPairing<DIM> pairings( dist_sources, target_centers(),
                                              d_use_knn, d_knn, d_radius );

    // Build the basis.
    Teuchos::RCP<Basis> basis = BP::create();

    // Build the interpolation matrix.
    Teuchos::ArrayRCP<SupportId> children_per_parent =
        pairings.childrenPerParent();
    SupportId max_entries_per_row = *std::max_element(
        children_per_parent.begin(), children_per_parent.end() );
    d_coupling_matrix = Teuchos::rcp( new Tpetra::CrsMatrix<Scalar, LO, GO>(
        range_map, max_entries_per_row ) );
    Teuchos::ArrayView<const double> target_view;
    Teuchos::Array<GO> indices( max_entries_per_row );
    Teuchos::ArrayView<const double> values;
    Teuchos::ArrayView<const unsigned> pair_gids;
    int nn = 0;
    int local_num_tgt = target_support_ids.size();
    for ( int i = 0; i < local_num_tgt; ++i )
    {
        // If there is no support for this target center then do not build a
        // local basis.
        if ( 0 < pairings.childCenterIds( i ).size() )
        {
            // Get a view of this target center.
            target_view = target_centers( i * DIM, DIM );

            // Build the local interpolation problem.
            LocalMLSProblem<Basis, DIM> local_problem(
                target_view, pairings.childCenterIds( i ), dist_sources, *basis,
                pairings.parentSupportRadius( i ) );

            // Get MLS shape function values for this target point.
            values = local_problem.shapeFunction();
            nn = values.size();

            // Populate the interpolation matrix row.
            pair_gids = pairings.childCenterIds( i );
            for ( int j = 0; j < nn; ++j )
            {
                indices[j] = dist_source_support_ids[pair_gids[j]];
            }
            d_coupling_matrix->insertGlobalValues( target_support_ids[i],
                                                   indices( 0, nn ), values );
        }
    }
    d_coupling_matrix->fillComplete( domain_map, range_map );
    DTK_ENSURE( d_coupling_matrix->isFillComplete() );
}

//---------------------------------------------------------------------------//
// Apply the operator.
template <class Basis, int DIM>
void MovingLeastSquareReconstructionOperator<Basis, DIM>::applyImpl(
    const TpetraMultiVector &X, TpetraMultiVector &Y, Teuchos::ETransp mode,
    double alpha, double beta ) const
{
    d_coupling_matrix->apply( X, Y, mode, alpha, beta );
}

//---------------------------------------------------------------------------//
// Transpose apply option.
template <class Basis, int DIM>
bool MovingLeastSquareReconstructionOperator<Basis,
                                             DIM>::hasTransposeApplyImpl() const
{
    return true;
}

//---------------------------------------------------------------------------//
// Extract node coordinates and ids from an iterator.
template <class Basis, int DIM>
void MovingLeastSquareReconstructionOperator<Basis, DIM>::getNodeCoordsAndIds(
    const Teuchos::RCP<FunctionSpace> &space, const int entity_dim,
    Teuchos::ArrayRCP<double> &centers,
    Teuchos::ArrayRCP<GO> &support_ids ) const
{
    // Get an iterator over the local nodes.
    EntityIterator iterator;
    if ( Teuchos::nonnull( space->entitySet() ) )
    {
        LocalEntityPredicate local_predicate(
            space->entitySet()->communicator()->getRank() );
        PredicateFunction predicate = PredicateComposition::And(
            space->selectFunction(), local_predicate.getFunction() );
        iterator = space->entitySet()->entityIterator( entity_dim, predicate );
    }

    // Extract the coordinates and support ids of the nodes.
    int local_num_node = iterator.size();
    centers = Teuchos::ArrayRCP<double>( DIM * local_num_node );
    support_ids = Teuchos::ArrayRCP<GO>( local_num_node );
    Teuchos::Array<SupportId> node_supports;
    EntityIterator begin = iterator.begin();
    EntityIterator end = iterator.end();
    int entity_counter = 0;
    for ( EntityIterator entity = begin; entity != end;
          ++entity, ++entity_counter )
    {
        space->shapeFunction()->entitySupportIds( *entity, node_supports );
        DTK_CHECK( 1 == node_supports.size() );
        support_ids[entity_counter] = node_supports[0];
        space->localMap()->centroid( *entity,
                                     centers( DIM * entity_counter, DIM ) );
    }
}

//---------------------------------------------------------------------------//

} // end namespace DataTransferKit

//---------------------------------------------------------------------------//

#endif // end DTK_MOVINGLEASTSQUARERECONSTRUCTIONOPERATOR_IMPL_HPP

//---------------------------------------------------------------------------//
// end DTK_MovingLeastSquareReconstructionOperator_impl.hpp
//---------------------------------------------------------------------------//
