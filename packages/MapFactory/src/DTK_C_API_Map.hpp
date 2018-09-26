/****************************************************************************
 * Copyright (c) 2012-2018 by the DataTransferKit authors                   *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the DataTransferKit library. DataTransferKit is     *
 * distributed under a BSD 3-clause license. For the licensing terms see    *
 * the LICENSE file in the top-level directory.                             *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
/*!
 * \file
 * \brief C adapter to DTK Maps.
 */
#ifndef DTK_C_API_MAP_HPP
#define DTK_C_API_MAP_HPP

#include <DTK_C_API.h>
#include <DTK_C_API.hpp>
#include <DTK_MovingLeastSquaresOperator.hpp>
#include <DTK_NearestNeighborOperator.hpp>
#include <DTK_ParallelTraits.hpp>
#include <DTK_PointCloudOperator.hpp>
#include <DTK_UserApplication.hpp>

#include <Teuchos_DefaultMpiComm.hpp>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <mpi.h>

#include <memory>
#include <string>
#include <tuple>

namespace DataTransferKit
{
//---------------------------------------------------------------------------//
// Map interface base class. This allows us to hide the device template
// parameter in MapImpl when we do the casting in the C interface
// implementation. We won't know what device type we have in a template when
// calling the C interface and we don't need to - the apply function
// implementation details will dispatch to the correct device type based on
// the construction.
struct DTK_Map
{
    virtual ~DTK_Map() = default;

    virtual void apply( const std::string &source_field_name,
                        const std::string &target_field_name ) = 0;
};

//---------------------------------------------------------------------------//
template <class MapExecSpace, class SourceMemSpace, class TargetMemSpace>
struct DTK_MapImpl : public DTK_Map
{
    using map_device_type = typename MapExecSpace::device_type;

    DTK_MapImpl( MPI_Comm comm, DTK_UserApplicationHandle source,
                 DTK_UserApplicationHandle target,
                 boost::property_tree::ptree const &ptree )
        : _source( reinterpret_cast<DTK_Registry *>( source )->_registry )
        , _target( reinterpret_cast<DTK_Registry *>( target )->_registry )
    {
        // Create a Teuchos comm from the MPI comm.
        auto teuchos_comm = Teuchos::rcp( new Teuchos::MpiComm<int>( comm ) );

        // FOR NOW JUST CREATE A NEAREST NEIGHBOR OPERATOR FOR DEMONSTRATION
        // PURPOSES. THIS WILL BE REPLACED BY A PROPER FACTORY.

        // Get coordinates from the source and target.
        auto source_nodes = _source.getNodeList().coordinates;
        auto target_nodes = _target.getNodeList().coordinates;

        // For now things are layout left in the interface so copy to a
        // matching layout that is compatible with the operator.
        Kokkos::View<Coordinate **, map_device_type> source_nodes_copy(
            "src_nodes_copy", source_nodes.extent( 0 ),
            source_nodes.extent( 1 ) );
        Kokkos::deep_copy( source_nodes_copy, source_nodes );
        Kokkos::View<Coordinate **, map_device_type> target_nodes_copy(
            "tgt_nodes_copy", target_nodes.extent( 0 ),
            target_nodes.extent( 1 ) );
        Kokkos::deep_copy( target_nodes_copy, target_nodes );

        auto const which_map = ptree.get<std::string>( "Map Type" );
        if ( which_map == "Nearest Neighbor" || which_map == "NN" )
            _map = std::unique_ptr<NearestNeighborOperator<map_device_type>>(
                new NearestNeighborOperator<map_device_type>(
                    teuchos_comm, source_nodes_copy, target_nodes_copy ) );
        else if ( which_map == "Moving Least Squares" || which_map == "MLS" )
        {
            auto const order = ptree.get<std::string>( "Order", "Linear" );
            if ( order == "Linear" || order == "1" )
                _map = std::unique_ptr<MovingLeastSquaresOperator<
                    map_device_type, Wendland<0>,
                    MultivariatePolynomialBasis<Linear, 3>>>(
                    new MovingLeastSquaresOperator<
                        map_device_type, Wendland<0>,
                        MultivariatePolynomialBasis<Linear, 3>>(
                        teuchos_comm, source_nodes_copy, target_nodes_copy ) );
            else if ( order == "Quadratic" || order == "2" )
                _map = std::unique_ptr<MovingLeastSquaresOperator<
                    map_device_type, Wendland<0>,
                    MultivariatePolynomialBasis<Quadratic, 3>>>(
                    new MovingLeastSquaresOperator<
                        map_device_type, Wendland<0>,
                        MultivariatePolynomialBasis<Quadratic, 3>>(
                        teuchos_comm, source_nodes_copy, target_nodes_copy ) );
            else
                throw std::runtime_error(
                    "Invalid order \"" + order +
                    "\" for creating a moving least squares map" );
        }
        else
            throw std::runtime_error( "Invalid map type \"" + which_map +
                                      "\"" );
    }

    void apply( const std::string &source_field_name,
                const std::string &target_field_name ) override
    {
        // Create fields.
        auto source_field = _source.getField( source_field_name );
        auto target_field = _target.getField( target_field_name );

        // Pull the data from the source.
        _source.pullField( source_field_name, source_field );

        // Copy to a compatible layout. Nearest neighbor operator needs a
        // compatible layout and only 1 dimension
        int num_src = source_field.dofs.extent( 0 );
        Kokkos::View<double *, map_device_type> source_field_copy(
            "source_field_copy", num_src );
        Kokkos::deep_copy(
            source_field_copy,
            Kokkos::subview( source_field.dofs, Kokkos::ALL, 0 ) );
        int num_tgt = target_field.dofs.extent( 0 );
        Kokkos::View<double *, map_device_type> target_field_copy(
            "target_field_copy", num_tgt );

        // Apply the map.
        _map->apply( source_field_copy, target_field_copy );

        // Copy the transferred field back to the original target layout.
        Kokkos::deep_copy( Kokkos::subview( target_field.dofs, Kokkos::ALL, 0 ),
                           target_field_copy );

        // Push the data to the target.
        _target.pushField( target_field_name, target_field );
    }

    UserApplication<double, SourceMemSpace> _source;
    UserApplication<double, TargetMemSpace> _target;
    std::unique_ptr<PointCloudOperator<map_device_type>> _map;
};

//---------------------------------------------------------------------------//
// Execution space validation.
bool validExecutionSpace( DTK_ExecutionSpace space )
{
    switch ( space )
    {
    case DTK_SERIAL:
#if defined( KOKKOS_ENABLE_SERIAL )
        return true;
        break;
#else
        return false;
        break;
#endif

    case DTK_OPENMP:
#if defined( KOKKOS_ENABLE_OPENMP )
        return true;
        break;
#else
        return false;
        break;
#endif

    case DTK_CUDA:
#if defined( KOKKOS_ENABLE_CUDA )
        return true;
        break;
#else
        return false;
        break;
#endif
    default:
        return false;
        break;
    }
}

// Memory space validation.
bool validMemorySpace( DTK_MemorySpace space )
{
    switch ( space )
    {
    case DTK_HOST_SPACE:
#if defined( KOKKOS_ENABLE_SERIAL ) || defined( KOKKOS_ENABLE_OPENMP )
        return true;
        break;
#else
        return false;
        break;
#endif
    case DTK_CUDAUVM_SPACE:
#if defined( KOKKOS_ENABLE_CUDA )
        return true;
        break;
#else
        return false;
        break;
#endif
    default:
        return false;
        break;
    }
}

// Map space validation.
bool validMapSpaces( DTK_ExecutionSpace map_space, DTK_MemorySpace source_space,
                     DTK_MemorySpace target_space )
{
    return validExecutionSpace( map_space ) &&
           validMemorySpace( source_space ) && validMemorySpace( target_space );
}

//---------------------------------------------------------------------------//
// Create a map.
DTK_Map *createMap( DTK_ExecutionSpace map_space, MPI_Comm comm,
                    DTK_UserApplicationHandle source,
                    DTK_UserApplicationHandle target, const char *options )
{
    // Parse options.
    std::stringstream ss;
    ss.str( options );
    boost::property_tree::ptree ptree;
    boost::property_tree::read_json( ss, ptree );

    // Get the map type.
    const std::string map_type( ptree.get<std::string>( "Map Type" ) );

    // Until a factory is added only nearest neighbor is supported.
    DTK_INSIST( "Nearest Neighbor" == map_type );

    // Get the user source and target memory spaces.
    DTK_MemorySpace src_space =
        reinterpret_cast<DataTransferKit::DTK_Registry *>( source )->_space;
    DTK_MemorySpace tgt_space =
        reinterpret_cast<DataTransferKit::DTK_Registry *>( target )->_space;

    // Check up front that we have been asked for execution and memory spaces
    // that are available in the kokkos build. This lets use a little cleaner
    // macro logic below because we know we will get a valid map space
    // already.
    DTK_INSIST( validMapSpaces( map_space, src_space, tgt_space ) );

    DTK_Map *map;
    switch ( map_space )
    {
    case DTK_SERIAL:
#if defined( KOKKOS_ENABLE_SERIAL )
        switch ( src_space )
        {
        case DTK_HOST_SPACE:
            switch ( tgt_space )
            {
            case DTK_HOST_SPACE:
                map = new DTK_MapImpl<Serial, HostSpace, HostSpace>(
                    comm, source, target, ptree );
                break;

            case DTK_CUDAUVM_SPACE:
#if defined( KOKKOS_ENABLE_CUDA )
                map = new DTK_MapImpl<Serial, HostSpace, CudaUVMSpace>(
                    comm, source, target, ptree );
#endif
                break;
            }
            break;

        case DTK_CUDAUVM_SPACE:
#if defined( KOKKOS_ENABLE_CUDA )
            switch ( tgt_space )
            {
            case DTK_HOST_SPACE:
                map = new DTK_MapImpl<Serial, CudaUVMSpace, HostSpace>(
                    comm, source, target, ptree );
                break;

            case DTK_CUDAUVM_SPACE:
                map = new DTK_MapImpl<Serial, CudaUVMSpace, CudaUVMSpace>(
                    comm, source, target, ptree );
                break;
            }
#endif
            break;
        }
#endif
        break;

    case DTK_OPENMP:
#if defined( KOKKOS_ENABLE_OPENMP )
        switch ( src_space )
        {
        case DTK_HOST_SPACE:
            switch ( tgt_space )
            {
            case DTK_HOST_SPACE:
                map = new DTK_MapImpl<OpenMP, HostSpace, HostSpace>(
                    comm, source, target, ptree );
                break;

            case DTK_CUDAUVM_SPACE:
#if defined( KOKKOS_ENABLE_CUDA )
                map = new DTK_MapImpl<OpenMP, HostSpace, CudaUVMSpace>(
                    comm, source, target, ptree );
#endif
                break;
            }
            break;

        case DTK_CUDAUVM_SPACE:
#if defined( KOKKOS_ENABLE_CUDA )
            switch ( tgt_space )
            {
            case DTK_HOST_SPACE:
                map = new DTK_MapImpl<OpenMP, CudaUVMSpace, HostSpace>(
                    comm, source, target, ptree );
                break;

            case DTK_CUDAUVM_SPACE:
                map = new DTK_MapImpl<OpenMP, CudaUVMSpace, CudaUVMSpace>(
                    comm, source, target, ptree );
                break;
            }
#endif
            break;
        }
#endif
        break;

    case DTK_CUDA:
#if defined( KOKKOS_ENABLE_CUDA )
        switch ( src_space )
        {
        case DTK_HOST_SPACE:
#if defined( KOKKOS_ENABLE_SERIAL ) || defined( KOKKOS_ENABLE_OPENMP )
            switch ( tgt_space )
            {
            case DTK_HOST_SPACE:
                map = new DTK_MapImpl<Cuda, HostSpace, HostSpace>(
                    comm, source, target, ptree );
                break;

            case DTK_CUDAUVM_SPACE:
                map = new DTK_MapImpl<Cuda, HostSpace, CudaUVMSpace>(
                    comm, source, target, ptree );
                break;
            }
#endif
            break;

        case DTK_CUDAUVM_SPACE:
            switch ( tgt_space )
            {
            case DTK_HOST_SPACE:
#if defined( KOKKOS_ENABLE_SERIAL ) || defined( KOKKOS_ENABLE_OPENMP )
                map = new DTK_MapImpl<Cuda, CudaUVMSpace, HostSpace>(
                    comm, source, target, ptree );
#endif
                break;
            case DTK_CUDAUVM_SPACE:
                map = new DTK_MapImpl<Cuda, CudaUVMSpace, CudaUVMSpace>(
                    comm, source, target, ptree );
                break;
            }
            break;
        }
#endif
        break;
    }

    return map;
}

//---------------------------------------------------------------------------//

} // namespace DataTransferKit

#endif // DTK_C_API_MAP_HPP
