/****************************************************************************
 * Copyright (c) 2012-2017 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 ****************************************************************************/

#ifndef DTK_INTERPOLATION_DECL_HPP
#define DTK_INTERPOLATION_DECL_HPP

#include "DTK_ConfigDefs.hpp"
#include <DTK_CellTypes.h>
#include <DTK_DetailsBox.hpp>
#include <DTK_DetailsPoint.hpp>
#include <DTK_DistributedSearchTree.hpp>
#include <DTK_FETypes.h>
#include <DTK_InterpolationOperator.hpp>
#include <DTK_QuadratureTypes.h>

#include <Shards_CellTopologyManagedData.hpp>
#include <Tpetra_Distributor.hpp>

#include <string>

namespace Teuchos
{
template <typename Ordinal>
class SerializationTraits<Ordinal, bool>
    : public DirectSerializationTraits<Ordinal, bool>
{
};
}

namespace DataTransferKit
{
template <typename DeviceType>
class Interpolation
{
  public:
    /**
     * Constructor. When this constructor is used getReferencePoints() can be
     * called but apply() cannot.
     * @param comm
     * @param cell_topologies_view
     * @param cells
     * @param nodes_coordinates coordinates of all the nodes in the mesh
     * @param points_coordinates coordinates in the physical frame of the points
     * that we want to evaluate.
     * @param fe_type finite element types used. The order is determined by
     * cell_topologies_view.
     */
    Interpolation(
        Teuchos::RCP<const Teuchos::Comm<int>> comm,
        Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view,
        Kokkos::View<unsigned int *, DeviceType> cells,
        Kokkos::View<double **, DeviceType> nodes_coordinates,
        Kokkos::View<double **, DeviceType> points_coordinates,
        DTK_FEType fe_type );

    /**
     * Constructor. When this constructor is used both getReferencePoints() and
     * apply() can be used.
     * @param comm
     * @param cell_topologies_view
     * @param cells
     * @param nodes_coordinates coordinates of all the nodes in the mesh
     * @param points_coordinates coordinates in the physical frame of the points
     * that we want to evaluate.
     * @param object_dof_ids dof ids associated to each node. This View is
     * similar to cell but instead of storing node ids, it stores dof ids.
     * @param fe_type finite element types used. The order is determined by
     * cell_topologies_view.
     */
    Interpolation(
        Teuchos::RCP<const Teuchos::Comm<int>> comm,
        Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view,
        Kokkos::View<unsigned int *, DeviceType> cells,
        Kokkos::View<double **, DeviceType> nodes_coordinates,
        Kokkos::View<double **, DeviceType> points_coordinates,
        Kokkos::View<LocalOrdinal *, DeviceType> object_dof_ids,
        DTK_FEType fe_type );

    /**
     * Constructor. When this constructor is used getReferencePoints() can be
     * called but apply() cannot.
     * @param comm
     * @param cell_topologies_view
     * @param cells
     * @param nodes_coordinates coordinates of all the nodes in the mesh
     * @param points_coordinates coordinates in the physical frame of the points
     * that we want to evaluate.
     * @param fe_order_view
     * @param quadrature_view
     * @param fe_type finite element types used. The order is determined by
     * cell_topologies_view.
     */
    Interpolation(
        Teuchos::RCP<const Teuchos::Comm<int>> comm,
        Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view,
        Kokkos::View<unsigned int *, DeviceType> cells,
        Kokkos::View<double **, DeviceType> nodes_coordinates,
        Kokkos::View<double **, DeviceType> points_coordinates,
        Kokkos::View<unsigned int *, DeviceType> fe_order_view,
        Kokkos::View<DTK_Quadrature *, DeviceType> quadrature_view,
        DTK_FEType fe_type );

    /**
     * Constructor. When this constructor is used getReferencePoints() can be
     * called but apply() cannot.
     * @param comm
     * @param cell_topologies_view
     * @param cells
     * @param nodes_coordinates coordinates of all the nodes in the mesh
     * @param points_coordinates coordinates in the physical frame of the points
     * that we want to evaluate.
     * @param object_dof_ids dof ids associated to each node. This View is
     * similar to cell but instead of storing node ids, it stores dof ids.
     * @param fe_order_view
     * @param quadrature_view
     * @param fe_type finite element types used. The order is determined by
     * cell_topologies_view.
     */
    Interpolation(
        Teuchos::RCP<const Teuchos::Comm<int>> comm,
        Kokkos::View<DTK_CellTopology *, DeviceType> cell_topologies_view,
        Kokkos::View<unsigned int *, DeviceType> cells,
        Kokkos::View<double **, DeviceType> nodes_coordinates,
        Kokkos::View<double **, DeviceType> points_coordinates,
        Kokkos::View<LocalOrdinal *, DeviceType> object_dof_ids,
        Kokkos::View<unsigned int *, DeviceType> fe_order_view,
        Kokkos::View<DTK_Quadrature *, DeviceType> quadrature_view,
        DTK_FEType fe_type );

    /**
     * This function performs the interpolation. It can be called only if
     * object_dofs_ids was given to the constructor.
     * @param [in] X (n dofs, n fields)
     * @param [out] Y (n phys points, n fields)
     */
    template <typename Scalar>
    void apply( Kokkos::View<Scalar **, DeviceType> X,
                Kokkos::View<Scalar **, DeviceType> Y );

    /**
     * This function gets:
     * @param [out] phys_points the coordinates of the points in the physical
     * frame (the points are not ordered)
     * @param [out] ref_points the coordinates of the points in the reference
     * frame
     * @param [out] pt_in_cell a flag to know if the point is in the cell
     *
     * All the Views will be resized by this function.
     */
    void getReferencePoints(
        Kokkos::View<DataTransferKit::Point *, DeviceType> &phys_points,
        Kokkos::View<DataTransferKit::Point *, DeviceType> &ref_points,
        Kokkos::View<bool *, DeviceType> &pt_in_cell );

    /**
     * Convert the 1D Kokkos View cells and coordinates to array of 3D Kokkos
     * View more suitable for Intrepid2 and the bounding boxes of each cell.
     *
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    void convertMesh( Kokkos::View<unsigned int *, DeviceType> cells,
                      Kokkos::View<double **, DeviceType> coordinates );

    /**
     * Perform the distributed search and sends the points and the cell indices
     * to the processors owning the cells.
     *
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    void performDistributedSearch(
        Kokkos::View<double **, DeviceType> points_coord,
        Kokkos::View<DataTransferKit::Point *, DeviceType> &imported_points,
        Kokkos::View<int *, DeviceType> &imported_query_ids,
        Kokkos::View<int *, DeviceType> &imported_cell_indices );

    /**
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    template <typename T1, typename T2>
    void computeOffset( Kokkos::View<T1 *, DeviceType> predicate, T2 value,
                        Kokkos::View<unsigned int *, DeviceType> offset );

    /**
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    void
    filterTopology( Kokkos::View<unsigned int *, DeviceType> topo,
                    unsigned int topo_id,
                    Kokkos::View<int *, DeviceType> cell_indices,
                    Kokkos::View<DataTransferKit::Point *, DeviceType> points,
                    Kokkos::View<int *, DeviceType> query_ids,
                    Kokkos::View<unsigned int *, DeviceType> map,
                    Kokkos::View<int *, DeviceType> filtered_cell_indices,
                    Kokkos::View<double **, DeviceType> filtered_points,
                    Kokkos::View<int *, DeviceType> filtered_query_ids );

    /**
     * @note This function should be <b>private</b> but lambda functions can
     * only be called from a public function in CUDA.
     */
    void findReferencePoints(
        Kokkos::View<unsigned int *, DeviceType> cells,
        Kokkos::View<double **, DeviceType> nodes_coordinates,
        Kokkos::View<double **, DeviceType> points_coordinates );

    /**
     * @note This function exists only because we cannot use lambda function in
     * the constructor with CUDA.
     */
    void computeNodeOffset(
        unsigned int const n_cells,
        Kokkos::View<DTK_CellTopology *, DeviceType> cell_topo_view,
        Kokkos::View<unsigned int[DTK_N_TOPO], DeviceType> n_nodes_per_topo,
        Kokkos::View<unsigned int *, DeviceType> nodes_per_cell,
        Kokkos::View<unsigned int *, DeviceType> node_offset );

    /**
     * @note This function exists only because we cannot lambda function in the
     * constructor with CUDA.
     */
    void
    fillCellDofIds( unsigned int const topo_id, unsigned int const n_cells,
                    unsigned int const n_nodes_per_topo,
                    Kokkos::View<DTK_CellTopology *, DeviceType> cell_topo_view,
                    Kokkos::View<unsigned int *, DeviceType> offset,
                    Kokkos::View<unsigned int *, DeviceType> node_offset,
                    Kokkos::View<LocalOrdinal *, DeviceType> object_dof_ids,
                    Kokkos::View<LocalOrdinal **, DeviceType> cell_dofs_ids );

  private:
    /**
     * Compute the number of cells of each topologies supported by DTK.
     */
    void computeTopologies();

    /**
     * Compute the coordinates in the reference frame and determine if the
     * points are in the cells.
     */
    void performPointInCell(
        Kokkos::View<double ***, DeviceType> cells,
        shards::CellTopology const &cell_topology,
        Kokkos::View<int *, DeviceType> imported_cell_indices,
        Kokkos::View<DataTransferKit::Point *, DeviceType> imported_points,
        Kokkos::View<int *, DeviceType> imported_query_ids,
        Kokkos::View<unsigned int *, DeviceType> topo, unsigned int topo_id,
        Kokkos::View<double **, DeviceType> filtered_points,
        Kokkos::View<int *, DeviceType> filtered_cell_indices,
        Kokkos::View<int *, DeviceType> filtered_query_ids,
        Kokkos::View<unsigned int *, DeviceType> points_indices_map,
        Kokkos::View<double **, DeviceType> reference_points,
        Kokkos::View<bool *, DeviceType> point_in_cell );

    /**
     * If true apply() can be called otherwise only getReferencePoints() can be
     * called.
     */
    bool _apply_allowed;
    /**
     * Dimension of domain (2 or 3).
     */
    unsigned int _dim;
    /**
     * Communicator.
     */
    Teuchos::RCP<const Teuchos::Comm<int>> _comm;
    /**
     * Distributor used to send the results back to the processor that started
     * the query.
     */
    Tpetra::Distributor _distributor;
    /**
     * Map between the topology enum and the shards objects.
     */
    std::array<shards::CellTopology, DTK_N_TOPO> _cell_topologies;
    /**
     * Number of cells of each of the supported topology present in the mesh.
     */
    std::array<unsigned int, DTK_N_TOPO> _n_cells_per_topo;
    /**
     * View of the ids associated to each query.
     */
    Kokkos::View<int *, DeviceType> _query_id;
    /**
     * Topology associated to each cell.
     */
    Kokkos::View<DTK_CellTopology *, DeviceType> _cell_topologies_view;
    /**
     * Array of Views of coordinates of the nodes of the cells. Each View in the
     * array is associated to one topology.
     */
    std::array<Kokkos::View<double ***, DeviceType>, DTK_N_TOPO> _block_cells;
    /**
     * Bounding boxes associated to the cells.
     */
    Kokkos::View<DataTransferKit::Box *, DeviceType> _bounding_boxes;
    /**
     * Map between bounding boxes and cells.
     */
    Kokkos::View<unsigned int **, DeviceType> _bounding_box_to_cell;
    /**
     * Map between the indices of the unfiltered points and the indices of the
     * filtered points.
     */
    std::array<Kokkos::View<unsigned int *, DeviceType>, DTK_N_TOPO>
        _point_indices_map;
    /**
     * Coordinates of the reference points (output of PointInCell).
     */
    std::array<Kokkos::View<double **, DeviceType>, DTK_N_TOPO>
        _reference_points;
    /**
     * Flags used to determine if the points are in the associated cells
     * (output of PointInCell).
     */
    std::array<Kokkos::View<bool *, DeviceType>, DTK_N_TOPO> _point_in_cell;
    /**
     * Points that are in the cells.
     */
    std::array<Kokkos::View<double **, DeviceType>, DTK_N_TOPO>
        _filtered_points;
    /**
     * Cells indices associated to _filtered_points.
     */
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO>
        _filtered_cell_indices;
    /**
     * Query ids associated to _filtered_points.
     */
    std::array<Kokkos::View<int *, DeviceType>, DTK_N_TOPO> _filtered_query_ids;
    /**
     * Dofs ids associated to each node of the cells.
     */
    std::array<Kokkos::View<LocalOrdinal **, DeviceType>, DTK_N_TOPO>
        _cell_dofs_ids;
};

template <typename DeviceType>
template <typename Scalar>
void Interpolation<DeviceType>::apply( Kokkos::View<Scalar **, DeviceType> X,
                                       Kokkos::View<Scalar **, DeviceType> Y )
{
    DTK_REQUIRE( _apply_allowed );
    DTK_REQUIRE( X.extent( 1 ) == Y.extent( 1 ) );
    // TODO what happens when a point is found in two cells? Or is on edge we
    // probably want to do sth different if we are using DG or continuous FEM
    using ExecutionSpace = typename DeviceType::execution_space;
    unsigned int dim = _dim;

    unsigned int n_exported_points = 0;
    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
        n_exported_points += _reference_points[topo_id].extent( 0 );

    unsigned int const n_fields = X.extent( 1 );
    Kokkos::View<Scalar **, DeviceType> exported_values(
        "exported_values", n_exported_points, n_fields );
    Kokkos::View<int *, DeviceType> exported_query_ids( "exported_query_ids",
                                                        n_exported_points );
    Kokkos::View<bool *, DeviceType> exported_pt_in_cell( "exported_pt_in_cell",
                                                          n_exported_points );

    for ( unsigned int topo_id = 0; topo_id < DTK_N_TOPO; ++topo_id )
    {
        unsigned int n_ref_points = _point_in_cell[topo_id].extent( 0 );
        if ( n_ref_points != 0 )
        {
            int n_filtered_ref_points = 0;
            Kokkos::View<bool *, DeviceType> pt_in_cell =
                _point_in_cell[topo_id];
            Kokkos::parallel_reduce(
                REGION_NAME( "compute_n_ref_pts" ),
                Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_points ),
                KOKKOS_LAMBDA( int i, int &partial_sum ) {
                    if ( pt_in_cell[i] == true )
                        partial_sum += 1;
                },
                n_filtered_ref_points );

            // We are only interested in points that belong to the cells. So we
            // need to filter out, all the points that were false positive of
            // the distributed search.
            Kokkos::View<Coordinate **, DeviceType> reference_points =
                _reference_points[topo_id];
            Kokkos::View<LocalOrdinal **, DeviceType> cell_dofs_ids =
                _cell_dofs_ids[topo_id];
            unsigned int const n_dofs_per_cell = cell_dofs_ids.extent( 1 );
            Kokkos::View<unsigned int *, DeviceType> filtered_map(
                "filtered_map", n_ref_points );
            Kokkos::View<Coordinate **, DeviceType> filtered_ref_points(
                "filtered_reference_points", n_filtered_ref_points, _dim );
            Kokkos::View<LocalOrdinal **, DeviceType> filtered_cell_dofs_ids(
                "filtered_cell_dofs_ids", n_filtered_ref_points,
                n_dofs_per_cell );
            Kokkos::View<int *, DeviceType> filtered_cell_indices =
                _filtered_cell_indices[topo_id];
            Kokkos::View<unsigned int *, DeviceType> offset( "offset",
                                                             n_ref_points );
            computeOffset( pt_in_cell, true, offset );
            Kokkos::parallel_for(
                REGION_NAME( "filter" ),
                Kokkos::RangePolicy<ExecutionSpace>( 0, n_ref_points ),
                KOKKOS_LAMBDA( int const i ) {
                    if ( pt_in_cell[i] )
                    {
                        unsigned int k = offset( i );
                        for ( unsigned int d = 0; d < dim; ++d )
                        {
                            filtered_ref_points( k, d ) =
                                reference_points( i, d );
                            filtered_map( i ) = k;
                        }
                        for ( unsigned int dof = 0; dof < n_dofs_per_cell;
                              ++dof )
                            filtered_cell_dofs_ids( k, dof ) = cell_dofs_ids(
                                filtered_cell_indices( i ), dof );
                    }
                } );
            Kokkos::fence();

            // Now that only, the points in the cells are left, we can do the
            // interpolation using Intrepid2.
            Kokkos::View<Scalar **, DeviceType> Y_topo(
                "Y_topo_" + std::to_string( topo_id ), n_filtered_ref_points,
                n_fields );
            Kokkos::deep_copy( Y_topo, 0. );
            unsigned int cell_topo_key = _cell_topologies[topo_id].getKey();
            if ( cell_topo_key ==
                 shards::getCellTopologyData<shards::Hexahedron<8>>()->key )
            {
                Functor::InterpolationOperator<
                    Scalar,
                    Intrepid2::Impl::Basis_HGRAD_HEX_C1_FEM::Serial<
                        Intrepid2::OPERATOR_VALUE>,
                    8, DeviceType>
                    interpolation_operator( filtered_ref_points,
                                            filtered_cell_dofs_ids, X, Y_topo );
                Kokkos::parallel_for( REGION_NAME( "interpolate_hex_8" ),
                                      Kokkos::RangePolicy<ExecutionSpace>(
                                          0, n_filtered_ref_points ),
                                      interpolation_operator );
            }
            else if ( cell_topo_key ==
                      shards::getCellTopologyData<shards::Hexahedron<27>>()
                          ->key )
            {
                Functor::InterpolationOperator<
                    Scalar,
                    Intrepid2::Impl::Basis_HGRAD_HEX_C2_FEM::Serial<
                        Intrepid2::OPERATOR_VALUE>,
                    27, DeviceType>
                    interpolation_operator( filtered_ref_points,
                                            filtered_cell_dofs_ids, X, Y_topo );
                Kokkos::parallel_for( REGION_NAME( "interpolate_hex_27" ),
                                      Kokkos::RangePolicy<ExecutionSpace>(
                                          0, n_filtered_ref_points ),
                                      interpolation_operator );
            }
            else if ( cell_topo_key ==
                      shards::getCellTopologyData<shards::Pyramid<5>>()->key )
            {
                Functor::InterpolationOperator<
                    Scalar,
                    Intrepid2::Impl::Basis_HGRAD_PYR_C1_FEM::Serial<
                        Intrepid2::OPERATOR_VALUE>,
                    5, DeviceType>
                    interpolation_operator( filtered_ref_points,
                                            filtered_cell_dofs_ids, X, Y_topo );
                Kokkos::parallel_for( REGION_NAME( "interpolate_pyr_5" ),
                                      Kokkos::RangePolicy<ExecutionSpace>(
                                          0, n_filtered_ref_points ),
                                      interpolation_operator );
            }
            else if ( cell_topo_key ==
                      shards::getCellTopologyData<shards::Quadrilateral<4>>()
                          ->key )
            {
                Functor::InterpolationOperator<
                    Scalar,
                    Intrepid2::Impl::Basis_HGRAD_QUAD_C1_FEM::Serial<
                        Intrepid2::OPERATOR_VALUE>,
                    4, DeviceType>
                    interpolation_operator( filtered_ref_points,
                                            filtered_cell_dofs_ids, X, Y_topo );
                Kokkos::parallel_for( REGION_NAME( "interpolate_quad_4" ),
                                      Kokkos::RangePolicy<ExecutionSpace>(
                                          0, n_filtered_ref_points ),
                                      interpolation_operator );
            }
            else if ( cell_topo_key ==
                      shards::getCellTopologyData<shards::Quadrilateral<9>>()
                          ->key )
            {
                Functor::InterpolationOperator<
                    Scalar,
                    Intrepid2::Impl::Basis_HGRAD_QUAD_C2_FEM::Serial<
                        Intrepid2::OPERATOR_VALUE>,
                    9, DeviceType>
                    interpolation_operator( filtered_ref_points,
                                            filtered_cell_dofs_ids, X, Y_topo );
                Kokkos::parallel_for( REGION_NAME( "interpolate_quad_9" ),
                                      Kokkos::RangePolicy<ExecutionSpace>(
                                          0, n_filtered_ref_points ),
                                      interpolation_operator );
            }
            else if ( cell_topo_key ==
                      shards::getCellTopologyData<shards::Tetrahedron<4>>()
                          ->key )
            {
                Functor::InterpolationOperator<
                    Scalar,
                    Intrepid2::Impl::Basis_HGRAD_TET_C1_FEM::Serial<
                        Intrepid2::OPERATOR_VALUE>,
                    4, DeviceType>
                    interpolation_operator( filtered_ref_points,
                                            filtered_cell_dofs_ids, X, Y_topo );
                Kokkos::parallel_for( REGION_NAME( "interpolate_tet_4" ),
                                      Kokkos::RangePolicy<ExecutionSpace>(
                                          0, n_filtered_ref_points ),
                                      interpolation_operator );
            }
            else if ( cell_topo_key ==
                      shards::getCellTopologyData<shards::Tetrahedron<10>>()
                          ->key )
            {
                Functor::InterpolationOperator<
                    Scalar,
                    Intrepid2::Impl::Basis_HGRAD_TET_C2_FEM::Serial<
                        Intrepid2::OPERATOR_VALUE>,
                    10, DeviceType>
                    interpolation_operator( filtered_ref_points,
                                            filtered_cell_dofs_ids, X, Y_topo );
                Kokkos::parallel_for( REGION_NAME( "interpolate_tet_10" ),
                                      Kokkos::RangePolicy<ExecutionSpace>(
                                          0, n_filtered_ref_points ),
                                      interpolation_operator );
            }
            else if ( cell_topo_key ==
                      shards::getCellTopologyData<shards::Triangle<3>>()->key )
            {
                Functor::InterpolationOperator<
                    Scalar,
                    Intrepid2::Impl::Basis_HGRAD_TRI_C1_FEM::Serial<
                        Intrepid2::OPERATOR_VALUE>,
                    3, DeviceType>
                    interpolation_operator( filtered_ref_points,
                                            filtered_cell_dofs_ids, X, Y_topo );
                Kokkos::parallel_for( REGION_NAME( "interpolate_tri_3" ),
                                      Kokkos::RangePolicy<ExecutionSpace>(
                                          0, n_filtered_ref_points ),
                                      interpolation_operator );
            }
            else if ( cell_topo_key ==
                      shards::getCellTopologyData<shards::Triangle<6>>()->key )
            {
                Functor::InterpolationOperator<
                    Scalar,
                    Intrepid2::Impl::Basis_HGRAD_TRI_C2_FEM::Serial<
                        Intrepid2::OPERATOR_VALUE>,
                    6, DeviceType>
                    interpolation_operator( filtered_ref_points,
                                            filtered_cell_dofs_ids, X, Y_topo );
                Kokkos::parallel_for( REGION_NAME( "interpolate_tri_6" ),
                                      Kokkos::RangePolicy<ExecutionSpace>(
                                          0, n_filtered_ref_points ),
                                      interpolation_operator );
            }
            else if ( cell_topo_key ==
                      shards::getCellTopologyData<shards::Wedge<6>>()->key )
            {
                Functor::InterpolationOperator<
                    Scalar,
                    Intrepid2::Impl::Basis_HGRAD_WEDGE_C1_FEM::Serial<
                        Intrepid2::OPERATOR_VALUE>,
                    6, DeviceType>
                    interpolation_operator( filtered_ref_points,
                                            filtered_cell_dofs_ids, X, Y_topo );
                Kokkos::parallel_for( REGION_NAME( "interpolate_wedge_6" ),
                                      Kokkos::RangePolicy<ExecutionSpace>(
                                          0, n_filtered_ref_points ),
                                      interpolation_operator );
            }
            else if ( cell_topo_key ==
                      shards::getCellTopologyData<shards::Wedge<18>>()->key )
            {
                Functor::InterpolationOperator<
                    Scalar,
                    Intrepid2::Impl::Basis_HGRAD_WEDGE_C2_FEM::Serial<
                        Intrepid2::OPERATOR_VALUE>,
                    18, DeviceType>
                    interpolation_operator( filtered_ref_points,
                                            filtered_cell_dofs_ids, X, Y_topo );
                Kokkos::parallel_for( REGION_NAME( "interpolate_wedge_18" ),
                                      Kokkos::RangePolicy<ExecutionSpace>(
                                          0, n_filtered_ref_points ),
                                      interpolation_operator );
            }
            else
            {
                throw std::runtime_error( "Not implemented" );
            }
            Kokkos::fence();

            // Because we want to reuse the distributor, we need to put the
            // results back in right order.
            Kokkos::View<unsigned int *, DeviceType> pt_indices_map =
                _point_indices_map[topo_id];

            Kokkos::View<int *, DeviceType> query_ids =
                _filtered_query_ids[topo_id];
            Kokkos::View<bool *, DeviceType> point_in_cell =
                _point_in_cell[topo_id];
            using ExecutionSpace = typename DeviceType::execution_space;
            unsigned int const n_topo_points = reference_points.extent( 0 );
            Kokkos::parallel_for(
                REGION_NAME( "merge_block_view" ),
                Kokkos::RangePolicy<ExecutionSpace>( 0, n_topo_points ),
                KOKKOS_LAMBDA( int const i ) {
                    for ( unsigned int j = 0; j < n_fields; ++j )
                        exported_values( pt_indices_map( i ), j ) =
                            Y_topo( filtered_map( i ), j );
                    exported_query_ids( pt_indices_map( i ) ) = query_ids( i );
                    exported_pt_in_cell( pt_indices_map( i ) ) =
                        point_in_cell( i );
                } );
            Kokkos::fence();
        }
    }

    // Communicate the results back to the calling processors
    unsigned int const n_imports = _distributor.getTotalReceiveLength();
    Kokkos::View<int *, DeviceType> imported_query_ids( "imported_query_ids",
                                                        n_imports );
    Kokkos::View<Scalar **, DeviceType> imported_values( "imported_values",
                                                         n_imports, n_fields );
    Kokkos::View<bool *, DeviceType> imported_pt_in_cell( "imported_pt_in_cell",
                                                          n_imports );

    DataTransferKit::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        _distributor, exported_query_ids, imported_query_ids );
    DataTransferKit::DistributedSearchTreeImpl<DeviceType>::sendAcrossNetwork(
        _distributor, exported_pt_in_cell, imported_pt_in_cell );

    // TODO change this: either create a new distributor or a new structure to
    // pack the values
    for ( unsigned int i = 0; i < n_fields; ++i )
    {
        Kokkos::View<Scalar *, DeviceType> exported_val( "exported_val",
                                                         n_exported_points );
        Kokkos::View<Scalar *, DeviceType> imported_val( "imported_val",
                                                         n_imports );

        Kokkos::parallel_for(
            REGION_NAME( "extract_exported_values" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_exported_points ),
            KOKKOS_LAMBDA( int const j ) {
                exported_val( j ) = exported_values( j, i );
            } );
        Kokkos::fence();

        DataTransferKit::DistributedSearchTreeImpl<
            DeviceType>::sendAcrossNetwork( _distributor, exported_val,
                                            imported_val );

        Kokkos::parallel_for(
            REGION_NAME( "copy_imported_values" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_imports ),
            KOKKOS_LAMBDA( int const j ) {
                imported_values( j, i ) = imported_val( j );
            } );
        Kokkos::fence();
    }

    // We know have the results on the processors that initiated the query but
    // we also have a lot of junk (points that are not in the cell) that we need
    // to filter.
    unsigned int n_in_cell_query_ids = 0;
    Kokkos::parallel_reduce(
        REGION_NAME( "compute_n_in_cell_query_ids" ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, n_imports ),
        KOKKOS_LAMBDA( int i, unsigned int &partial_sum ) {
            if ( imported_pt_in_cell[i] == true )
                partial_sum += 1;
        },
        n_in_cell_query_ids );

    Kokkos::View<Scalar **, DeviceType> in_cell_values(
        "in_cell_values", n_in_cell_query_ids, n_fields );
    Kokkos::View<int *, DeviceType> in_cell_query_ids( "in_cell_query_ids",
                                                       n_in_cell_query_ids );

    Kokkos::View<unsigned int *, DeviceType> in_cell_offset( "in_cell_offset",
                                                             n_imports );
    computeOffset( imported_pt_in_cell, true, in_cell_offset );
    Kokkos::parallel_for( REGION_NAME( "filter_in_cell" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, n_imports ),
                          KOKKOS_LAMBDA( int const i ) {
                              if ( imported_pt_in_cell[i] )
                              {
                                  unsigned int k = in_cell_offset( i );
                                  for ( unsigned int j = 0; j < n_fields; ++j )
                                      in_cell_values( k, j ) =
                                          imported_values( i, j );
                                  in_cell_query_ids( k ) =
                                      imported_query_ids( i );
                              }
                          } );
    Kokkos::fence();

    if ( n_in_cell_query_ids != 0 )
    {
        // Because of the MPI communications and the sorting by topologies, all
        // the query have been reordered. So we put them back in the initial
        // order using the query ids.
        typedef Kokkos::BinOp1D<Kokkos::View<int *, DeviceType>> CompType;

        Kokkos::Experimental::MinMaxScalar<int> result;
        Kokkos::Experimental::MinMax<int> reducer( result );
        parallel_reduce(
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_in_cell_query_ids ),
            Kokkos::Impl::min_max_functor<Kokkos::View<int *, DeviceType>>(
                in_cell_query_ids ),
            reducer );
        if ( result.min_val == result.max_val )
            return;
        Kokkos::BinSort<Kokkos::View<int *, DeviceType>, CompType> bin_sort(
            in_cell_query_ids,
            CompType( n_in_cell_query_ids / 2, result.min_val, result.max_val ),
            true );
        bin_sort.create_permute_vector();
        bin_sort.sort( in_cell_query_ids );
        bin_sort.sort( in_cell_values );

        Kokkos::View<unsigned int *, DeviceType> mask( "mask",
                                                       n_in_cell_query_ids );
        Kokkos::deep_copy( mask, 1 );
        Kokkos::parallel_for(
            REGION_NAME( "compute_mask" ),
            Kokkos::RangePolicy<ExecutionSpace>( 1, n_in_cell_query_ids ),
            KOKKOS_LAMBDA( int const i ) {
                if ( in_cell_query_ids( i - 1 ) == in_cell_query_ids( i ) )
                    mask( i ) = 0;
            } );
        Kokkos::fence();

        // We have finally all the values in the right order but before we can
        // finally return the values we need to filter the data one more time.
        // Some points are correctly found on multiple cells, e.g. point on
        // vertices, so we need to get ride of the duplicates.
        Kokkos::View<unsigned int *, DeviceType> query_offset(
            "query_offset", n_in_cell_query_ids );
        Kokkos::parallel_scan(
            REGION_NAME( "compute_query_offset" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_in_cell_query_ids ),
            KOKKOS_LAMBDA( int i, int &update, bool final_pass ) {
                if ( final_pass )
                    query_offset( i ) = update;
                update += mask( i );
            } );
        Kokkos::fence();

        Kokkos::parallel_for(
            REGION_NAME( "fill_Y(0,i)" ),
            Kokkos::RangePolicy<ExecutionSpace>( 0, n_fields ),
            KOKKOS_LAMBDA( int const i ) {
                Y( 0, i ) = in_cell_values( 0, i );
            } );
        Kokkos::fence();

        Kokkos::parallel_for(
            REGION_NAME( "fill_Y" ),
            Kokkos::RangePolicy<ExecutionSpace>( 1, n_in_cell_query_ids ),
            KOKKOS_LAMBDA( int const i ) {
                if ( in_cell_query_ids( i - 1 ) != in_cell_query_ids( i ) )
                {
                    unsigned int k = query_offset( i );
                    for ( unsigned int j = 0; j < n_fields; ++j )
                        Y( k, j ) = in_cell_values( i, j );
                }
            } );
        Kokkos::fence();
    }
}

template <typename DeviceType>
template <typename T1, typename T2>
void Interpolation<DeviceType>::computeOffset(
    Kokkos::View<T1 *, DeviceType> predicate, T2 value,
    Kokkos::View<unsigned int *, DeviceType> offset )
{
    // Create a Kokkos::View with ones where predicate matches value and
    // with zeros everywhere else.
    unsigned int const size = predicate.extent( 0 );
    using ExecutionSpace = typename DeviceType::execution_space;
    Kokkos::View<unsigned int *, DeviceType> mask( "mask", size );
    Kokkos::parallel_for( REGION_NAME( "compute_mask" ),
                          Kokkos::RangePolicy<ExecutionSpace>( 0, size ),
                          KOKKOS_LAMBDA( int const i ) {
                              if ( predicate( i ) == value )
                                  mask( i ) = 1;
                              else
                                  mask( i ) = 0;
                          } );
    Kokkos::fence();

    // Compute an offset that is used be fill filtered_points
    Kokkos::parallel_scan(
        REGION_NAME( "compute_offset" ),
        Kokkos::RangePolicy<ExecutionSpace>( 0, size ),
        KOKKOS_LAMBDA( int i, int &update, bool final_pass ) {
            if ( final_pass )
                offset( i ) = update;
            update += mask( i );
        } );
    Kokkos::fence();
}
}

#endif
