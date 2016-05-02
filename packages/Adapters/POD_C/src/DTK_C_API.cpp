#include <DTK_C_API.h>
#include <DTK_PointCloudOperatorFactory.hpp>
#include <DTK_Point.hpp>
#include <DTK_BasicGeometryManager.hpp>
#include <DTK_EntityCenteredField.hpp>
#include <DTK_FieldMultiVector.hpp>
#include <DTK_DBC.hpp>

#include <Teuchos_DefaultMpiComm.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <sstream> 
#include <iostream>
#include <string>

//----------------------------------------------------------------------------//
Teuchos::RCP<Tpetra::Map<int,DataTransferKit::EntityId> const>
  build_contiguous_map(
    Teuchos::RCP<Teuchos::Comm<int> const> const& comm,
    DataTransferKit::EntityId local_num )
{
  DataTransferKit::EntityId global_num = 0;
  Teuchos::reduceAll<int,DataTransferKit::EntityId>( *comm,
                                                     Teuchos::REDUCE_SUM,
                                                     1,
                                                     &local_num,
                                                     &global_num );
  return Tpetra::createContigMap<int, DataTransferKit::EntityId>(
    global_num, local_num, comm );
}

//----------------------------------------------------------------------------//
Teuchos::Array<DataTransferKit::Entity>
  pull_entities(
    Teuchos::RCP<Teuchos::Comm<int> const> const& comm,
    double const* points,
    unsigned local_num,
    DTK_Data_layout data_layout,
    int space_dim )
{
  Teuchos::Array<DataTransferKit::Entity> entities(local_num);
  Teuchos::Array<double> coord(space_dim);
  DataTransferKit::EntityId global_id;
  int const comm_size = comm->getSize();
  int const comm_rank = comm->getRank();

  // @Stuart: any way to do that with some teuchos rountine?
  // I guess the tpetra map can give this information
  // Then I'd pass the map instead of comm and local_num
  std::vector<unsigned> shift(comm_size);
  Teuchos::gatherAll<int,unsigned>( *comm,
                                    1,
                                    &local_num,
                                    shift.size(),
                                    shift.data() );
  for (std::size_t k = 1; k < comm_size; ++k)
    shift[k] = shift[k-1];
  shift[0] = 0;
  for (std::size_t k = 1; k < comm_size; ++k)
    shift[k] += shift[k-1];

  // TODO: might be worth implementing some data extractor to be used both here
  // and when handling fields
  if (data_layout == DTK_INTERLEAVED)
    for (unsigned i = 0; i < local_num; ++i)
    {
      global_id = shift[comm_rank] + i;
      for (int j = 0; j < space_dim; ++j)
        coord[j] = points[space_dim*i+j];
      entities[i] = DataTransferKit::Point(global_id, comm_rank, coord);
    }
  else if (data_layout == DTK_BLOCKED)
    for (unsigned i = 0; i < local_num; ++i)
    {
      global_id = shift[comm_rank] + i;
      for (int j = 0; j < space_dim; ++j)
        coord[j] = points[i+local_num*j];
      entities[i] = DataTransferKit::Point(global_id, comm_rank, coord);
    }
  else
    throw std::runtime_error("Invalid data layout");
  return entities;
}

//----------------------------------------------------------------------------//
DTK_Map* DTK_Map_create( MPI_Comm        comm,
                         double const*   src_coord,
                         unsigned        src_num,
                         DTK_Data_layout src_layout,
                         double const*   tgt_coord,
                         unsigned        tgt_num,
                         DTK_Data_layout tgt_layout,
                         int             space_dim,
                         char const*     options )
{
  // Parse the options and build the parameter list
  int comm_rank;
  MPI_Comm_rank(comm, &comm_rank);
  if (comm_rank == 0)
    std::cout << "\nDTK_Map_create called with options:\n" << options << "\n";
  std::stringstream ss;
  ss.str(options);
  boost::property_tree::ptree ptree;
  boost::property_tree::read_json(ss, ptree);
  Teuchos::ParameterList parameters;
  parameters.set("Map Type", ptree.get<std::string>("Map Type",
    "Moving Least Square Reconstruction") );
  parameters.set("Spatial Dimension", space_dim);
  parameters.set("Basis Type", ptree.get<std::string>("Basis Type",
    "Wendland"));
  parameters.set("Basis Order", ptree.get<int>("Basis Order", 2));
  // TODO: kNN
  parameters.set("RBF Radius", ptree.get<double>("RBF Radius"));

  // Create the domain and range maps
  Teuchos::RCP<Teuchos::Comm<int>> teuchos_comm =
    Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(comm));

  auto domain_map = build_contiguous_map(teuchos_comm, src_num);
  auto range_map = build_contiguous_map(teuchos_comm, tgt_num);

  // Build the actual DTK map operator
  auto factory = DataTransferKit::PointCloudOperatorFactory();
  auto map_operator = factory.create( domain_map,
                                      range_map,
                                      parameters );

  // Set it up
  auto domain_entities = pull_entities(
    teuchos_comm, src_coord, src_num, src_layout, space_dim );
  auto range_entities = pull_entities(
    teuchos_comm, tgt_coord, tgt_num, tgt_layout, space_dim );

  DataTransferKit::BasicGeometryManager domain_manager(
    teuchos_comm, space_dim, domain_entities() );
  DataTransferKit::BasicGeometryManager range_manager(
    teuchos_comm, space_dim, range_entities() );

  map_operator->setup( domain_manager.functionSpace(),
                       range_manager.functionSpace() );

  // Return an opaque pointer. User is responsible for calling delete_map(...)
  return static_cast<DTK_Map*>(map_operator.release().get());
}


//----------------------------------------------------------------------------//
void DTK_Map_apply( DTK_Map*        dtk_map,
                    double const*   src_data,
                    DTK_Data_layout src_layout,
                    double*         tgt_data,
                    DTK_Data_layout tgt_layout,
                    int             field_dim )
{
  // Cast the opaque pointer back to a DTK map operator
  auto map_operator = static_cast<DataTransferKit::MapOperator*>(dtk_map);

  // Helper function to map data layouts
  auto dtk_entity_centered_field_layout = [](DTK_Data_layout data_layout)
    {
      if (data_layout == DTK_BLOCKED)
        return DataTransferKit::EntityCenteredField::BLOCKED;
      else if (data_layout == DTK_INTERLEAVED)
        return DataTransferKit::EntityCenteredField::INTERLEAVED;
      else
        throw std::runtime_error("Invalid data layout");
    };

  // Wrap the source vector
  auto domain_map = map_operator->getDomainMap();

  DataTransferKit::EntityCenteredField src_field(
    domain_map->getNodeElementList(),
    field_dim,
    Teuchos::ArrayRCP<double>(
      const_cast<double*>(src_data),
      0,
      field_dim * domain_map->getNodeNumElements(),
      false ),
    dtk_entity_centered_field_layout(src_layout) );

  DataTransferKit::FieldMultiVector domain_vector(
    domain_map->getComm(), Teuchos::rcpFromRef(src_field) );

  // Wrap the target vector
  auto range_map = map_operator->getRangeMap();

  DataTransferKit::EntityCenteredField tgt_field(
    range_map->getNodeElementList(),
    field_dim,
    Teuchos::ArrayRCP<double>(
      tgt_data,
      0,
      field_dim * range_map->getNodeNumElements(),
      false ),
    dtk_entity_centered_field_layout(tgt_layout) );

  DataTransferKit::FieldMultiVector range_vector(
    range_map->getComm(), Teuchos::rcpFromRef(tgt_field) );

  // Apply the map operator
  map_operator->apply(domain_vector, range_vector);
}

//----------------------------------------------------------------------------//
void DTK_Map_delete( DTK_Map * dtk_map )
{
  delete static_cast<DataTransferKit::MapOperator*>(dtk_map);
}
