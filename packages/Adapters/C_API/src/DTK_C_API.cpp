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
#include <DTK_C_API.h>
#include <DTK_DBC.hpp>
#include <DTK_EntityCenteredField.hpp>
#include <DTK_EntityCenteredShapeFunction.hpp>
#include <DTK_FieldMultiVector.hpp>
#include <DTK_MapOperator.hpp>
#include <DTK_POD_PointCloudEntitySet.hpp>
#include <DTK_POD_PointCloudLocalMap.hpp>
#include <DTK_POD_Types.hpp>
#include <DTK_PointCloudOperatorFactory.hpp>

#include <Teuchos_DefaultMpiComm.hpp>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <sstream>
#include <string>

//----------------------------------------------------------------------------//
Teuchos::RCP<const typename DataTransferKit::MapOperator::TpetraMap>
build_contiguous_map( Teuchos::RCP<Teuchos::Comm<int> const> const &comm,
                      DataTransferKit::EntityId local_num )
{
    DataTransferKit::EntityId global_num = 0;
    Teuchos::reduceAll<int, DataTransferKit::EntityId>(
        *comm, Teuchos::REDUCE_SUM, 1, &local_num, &global_num );
    return Tpetra::createContigMapWithNode<
        typename DataTransferKit::MapOperator::LO,
        typename DataTransferKit::MapOperator::GO,
        typename DataTransferKit::MapOperator::Node>( global_num, local_num,
                                                      comm );
}

//---------------------------------------------------------------------------//
DataTransferKit::DataLayout getLayout( DTK_Data_layout layout )
{
    if ( DTK_BLOCKED == layout )
    {
        return DataTransferKit::BLOCKED;
    }
    else if ( DTK_INTERLEAVED == layout )
    {
        return DataTransferKit::INTERLEAVED;
    }
    else
    {
        // Will always throw if we get a bad layout.
        DTK_INSIST( DTK_BLOCKED == layout || DTK_INTERLEAVED == layout );
        return DataTransferKit::BLOCKED;
    }
}

//---------------------------------------------------------------------------//
Teuchos::RCP<DataTransferKit::FunctionSpace> createFunctionSpace(
    const Teuchos::RCP<const Teuchos::Comm<int>> &comm, const double *coord,
    const Teuchos::ArrayView<const DataTransferKit::EntityId> &global_ids,
    DTK_Data_layout layout, unsigned num, int space_dim )
{
    // Get the data layout
    DataTransferKit::DataLayout data_layout = getLayout( layout );

    // Build the entity set.
    auto entity_set =
        Teuchos::rcp( new DataTransferKit::POD_PointCloudEntitySet(
            comm, coord, global_ids.getRawPtr(), num, space_dim,
            data_layout ) );

    // Build the local map.
    auto local_map =
        Teuchos::rcp( new DataTransferKit::POD_PointCloudLocalMap() );

    // Build the shape function.
    auto shape_function =
        Teuchos::rcp( new DataTransferKit::EntityCenteredShapeFunction() );

    // Build the function space.
    return Teuchos::rcp( new DataTransferKit::FunctionSpace(
        entity_set, local_map, shape_function, Teuchos::null ) );
}

//----------------------------------------------------------------------------//
DTK_Map *DTK_Map_create_f( MPI_Fint fint, double const *src_coord,
                           unsigned src_num, DTK_Data_layout src_layout,
                           double const *tgt_coord, unsigned tgt_num,
                           DTK_Data_layout tgt_layout, int space_dim,
                           char const *options )
{
    MPI_Comm comm = MPI_Comm_f2c( fint );
    return DTK_Map_create( comm, src_coord, src_num, src_layout, tgt_coord,
                           tgt_num, tgt_layout, space_dim, options );
}

//----------------------------------------------------------------------------//
DTK_Map *DTK_Map_create( MPI_Comm comm, double const *src_coord,
                         unsigned src_num, DTK_Data_layout src_layout,
                         double const *tgt_coord, unsigned tgt_num,
                         DTK_Data_layout tgt_layout, int space_dim,
                         char const *options )
{
    // Parse the options and build the parameter list
    std::stringstream ss;
    ss.str( options );
    boost::property_tree::ptree ptree;
    boost::property_tree::read_json( ss, ptree );
    Teuchos::ParameterList parameters;
    std::string const map_type( ptree.get<std::string>(
        "Map Type", "Moving Least Square Reconstruction" ) );
    parameters.set( "Map Type", map_type );
    parameters.set( "Spatial Dimension", space_dim );
    if ( map_type.compare( "Node To Node" ) != 0 )
    {
        parameters.set( "Basis Type",
                        ptree.get<std::string>( "Basis Type", "Wendland" ) );
        parameters.set( "Basis Order", ptree.get<int>( "Basis Order", 2 ) );
        std::string const search_type =
            ptree.get<std::string>( "Search Type", "Radius" );
        parameters.set( "Type of Search", search_type );
        if ( search_type.compare( "Radius" ) == 0 )
        {
            parameters.set( "RBF Radius", ptree.get<double>( "RBF Radius" ) );
        }
        else if ( search_type.compare( "Nearest Neighbor" ) == 0 )
        {
            parameters.set( "Num Neighbors",
                            ptree.get<int>( "Num Neighbors" ) );
        }
    }

    // Wrap the communicator.
    Teuchos::RCP<const Teuchos::Comm<int>> teuchos_comm =
        Teuchos::createMpiComm<int>( Teuchos::opaqueWrapper<MPI_Comm>( comm ) );

    // Create the domain and range maps
    auto domain_map = build_contiguous_map( teuchos_comm, src_num );
    auto range_map = build_contiguous_map( teuchos_comm, tgt_num );

    // Build the actual DTK map operator
    auto factory = DataTransferKit::PointCloudOperatorFactory();
    auto map_operator = factory.create( domain_map, range_map, parameters );

    // Create the function spaces.
    auto domain_space = createFunctionSpace( teuchos_comm, src_coord,
                                             domain_map->getNodeElementList(),
                                             src_layout, src_num, space_dim );

    auto range_space = createFunctionSpace( teuchos_comm, tgt_coord,
                                            range_map->getNodeElementList(),
                                            tgt_layout, tgt_num, space_dim );

    // Set up the map.
    map_operator->setup( domain_space, range_space );

    // Return an opaque pointer. User is responsible for calling delete_map(...)
    return static_cast<DTK_Map *>( map_operator.release().get() );
}

//----------------------------------------------------------------------------//
void DTK_Map_apply( DTK_Map *dtk_map, double const *src_data,
                    DTK_Data_layout src_layout, double *tgt_data,
                    DTK_Data_layout tgt_layout, int field_dim, bool transpose )
{
    // Cast the opaque pointer back to a DTK map operator
    auto map_operator = static_cast<DataTransferKit::MapOperator *>( dtk_map );

    // Helper function to map data layouts
    auto dtk_entity_centered_field_layout = []( DTK_Data_layout data_layout ) {
        if ( data_layout == DTK_BLOCKED )
            return DataTransferKit::EntityCenteredField::BLOCKED;
        else if ( data_layout == DTK_INTERLEAVED )
            return DataTransferKit::EntityCenteredField::INTERLEAVED;
        else
            throw std::runtime_error( "Invalid data layout" );
    };

    // Wrap the source vector
    auto domain_map = map_operator->getDomainMap();

    DataTransferKit::EntityCenteredField src_field(
        domain_map->getNodeElementList(), field_dim,
        Teuchos::ArrayRCP<double>( const_cast<double *>( src_data ), 0,
                                   field_dim * domain_map->getNodeNumElements(),
                                   false ),
        dtk_entity_centered_field_layout( src_layout ) );

    DataTransferKit::FieldMultiVector domain_vector(
        domain_map->getComm(), Teuchos::rcpFromRef( src_field ) );

    // Wrap the target vector
    auto range_map = map_operator->getRangeMap();

    DataTransferKit::EntityCenteredField tgt_field(
        range_map->getNodeElementList(), field_dim,
        Teuchos::ArrayRCP<double>(
            tgt_data, 0, field_dim * range_map->getNodeNumElements(), false ),
        dtk_entity_centered_field_layout( tgt_layout ) );

    DataTransferKit::FieldMultiVector range_vector(
        range_map->getComm(), Teuchos::rcpFromRef( tgt_field ) );

    // Apply the map operator
    map_operator->apply( domain_vector, range_vector,
                         transpose ? Teuchos::TRANS : Teuchos::NO_TRANS );
}

//----------------------------------------------------------------------------//
void DTK_Map_delete( DTK_Map *dtk_map )
{
    delete static_cast<DataTransferKit::MapOperator *>( dtk_map );
}

//---------------------------------------------------------------------------//
// end DTK_C_API.cpp
//---------------------------------------------------------------------------//
