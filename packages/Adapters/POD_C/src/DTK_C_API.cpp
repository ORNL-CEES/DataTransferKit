#include <DTK_C_API.h>
#include <DTK_PointCloudOperatorFactory.hpp>
#include <DTK_Point.hpp>
#include <DTK_BasicGeometryManager.hpp>
#include <DTK_DBC.hpp>

#include <Teuchos_DefaultMpiComm.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

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
Teuchos::Array<double>
  pull_field(
    double const* data,
    unsigned local_num,
    DTK_Data_layout data_layout,
    int field_dim )
{
  Teuchos::Array<double> field(local_num*field_dim);
  if (data_layout == DTK_INTERLEAVED)
    for (unsigned i = 0; i < local_num; ++i)
      for (int j = 0; j < field_dim; ++j)
        field[i+local_num*j] = data[field_dim*i+j];
  else if (data_layout == DTK_BLOCKED)
    std::copy(data, data+field_dim*local_num, field.getRawPtr());
  else
    throw std::runtime_error("Invalid data layout");
  return field;
}

//----------------------------------------------------------------------------//
void push_field(
  Teuchos::Array<double> const& field,
  double * data,
  DTK_Data_layout data_layout,
  int field_dim )
{
  unsigned const local_num = field.size() / field_dim;
  if (data_layout == DTK_INTERLEAVED)
    for (unsigned i = 0; i < local_num; ++i)
      for (int j = 0; j < field_dim; ++j)
        data[field_dim*i+j] = field[i+local_num*j];
  else if (data_layout == DTK_BLOCKED)
    std::copy(field.begin(), field.end(), data);
  else
    throw std::runtime_error("Invalid data layout");
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
  boost::property_tree::read_info(ss, ptree);
  Teuchos::ParameterList parameters;
  parameters.set("Map Type", ptree.get<std::string>("Map Type",
    "Moving Least Square Reconstruction") );
  parameters.set("Spatial Dimension", space_dim);
  parameters.set("Basis Type", ptree.get<std::string>("Basis Type", "Wu"));
  parameters.set("Basis Order", ptree.get<int>("Basis Order", 0));
  // TODO: kNN
  parameters.set("RBF Radius", ptree.get<double>("RBF Radius"));

  // Create the domain and range maps
  // @Stuart: Not sure that's what I'm supposed to do to get the communicator
  // from the raw MPI comm
  Teuchos::RCP<Teuchos::Comm<int>> teuchos_comm =
    Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(comm));

  auto domain_map = build_contiguous_map(teuchos_comm, src_num);
  auto range_map = build_contiguous_map(teuchos_comm, tgt_num);

  // Build the actual DTK map operator
  // TODO: PointCloudOperatorFactory::create(...) should be a static method!
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
                    double const*   src_field,
                    DTK_Data_layout src_layout,
                    double*         tgt_field,
                    DTK_Data_layout tgt_layout,
                    int             field_dim )
{
  // Cast the opaque pointer back to a DTK map operator
  auto map_operator = static_cast<DataTransferKit::MapOperator*>(dtk_map);

  // Get a communicator
  auto domain_comm = map_operator->getDomainMap()->getComm();
  int const comm_rank = domain_comm->getRank();

  // Wrap the source and target vectors
  // TODO:
  // @Stuart: Can you help with that?
  auto domain_map = map_operator->getDomainMap();
  DTK_CHECK( domain_map->isContiguous() );
  unsigned src_local_num = domain_map->getMaxLocalIndex() - domain_map->getMinLocalIndex() + 1;
  Teuchos::Array<double> range_field = pull_field(
    src_field, src_local_num, src_layout, field_dim );
  DTK_CHECK( range_field.size() == src_local_num*field_dim );

  auto range_map = map_operator->getRangeMap();
  DTK_CHECK( range_map->isContiguous() );
  unsigned tgt_local_num = range_map->getMaxLocalIndex() - range_map->getMinLocalIndex() + 1;
  Teuchos::Array<double> domain_field(tgt_local_num*field_dim);

  // TODO: I suspect I need the domain and range managers to build the
  // vectors...
  // In that case I will have to define some kind of struct that keeps them
  // around next to the map operator :(

  //auto domain_vector = DataTransferKit::FieldMultiVector(
  //  Teuchos::rcp(new DataTransferKit::EntityCenteredField(domain_entities(), field_dim, domain_field)),
  //  domain_manager.functionSpace()->entitySet() );

  // Apply the map operator
  // TODO:
  // map_operator->apply(domain_vector, range_vector);

  push_field(domain_field, tgt_field, tgt_layout, field_dim);
}

//----------------------------------------------------------------------------//
void DTK_Map_delete( DTK_Map * dtk_map )
{
  delete static_cast<DataTransferKit::MapOperator*>(dtk_map);
}
