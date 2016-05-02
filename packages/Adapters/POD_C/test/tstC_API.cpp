#include "DTK_C_API.h"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "mpi.h"

namespace {

TEUCHOS_UNIT_TEST( C_API, create_apply_delete_map )
{
  // get the raw mpi communicator
  Teuchos::RCP<const Teuchos::Comm<int>> teuchos_comm =
    Teuchos::DefaultComm<int>::getComm();
  MPI_Comm mpi_comm =
    *Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(
      teuchos_comm )->getRawMpiComm();

  // Using JSON format
  std::string const options = "{ "\
    "\"Map Type\": \"Spline Interpolation\", "\
    "\"Basis Type\": \"Wendland\", "\
    "\"Basis Order\": 0, "\
    "\"RBF Radius\": 0.1 }";

  // TODO: fill these
  std::srand(0);
  int const space_dim = 3;
  int const field_dim = 1;
  unsigned const src_num = 100;
  std::vector<double> src_coord(space_dim*src_num);
  std::vector<double> src_field(field_dim*src_num);
  for (unsigned i; i < src_num; ++i)
  {
    for (unsigned j; j < space_dim; ++j)
    {
      // assume data layout is contiguous
      src_coord[i+j*src_num] = std::rand() / RAND_MAX;
    }
    for (unsigned j; j < field_dim; ++j)
    {
      // assume data layout is interleaved
      src_field[i*field_dim+j] = static_cast<double>(j);
    }
  }
  unsigned const tgt_num = 50;
  std::vector<double> tgt_coord(space_dim*tgt_num);
  std::vector<double> tgt_field(field_dim*tgt_num);
  for (unsigned i; i < tgt_num; ++i)
  {
    for (unsigned j; j < space_dim; ++j)
    {
      tgt_coord[i+j*tgt_num] = std::rand() / RAND_MAX;
    }
    for (unsigned j; j < field_dim; ++j)
    {
      tgt_field[i*field_dim+j] = static_cast<double>(j);
    }
  }
  TEST_EQUALITY( src_coord.size(), space_dim * src_num );
  TEST_EQUALITY( tgt_coord.size(), space_dim * tgt_num );

  DTK_Map* dtk_map = DTK_Map_create( mpi_comm,
                                     src_coord.data(),
                                     src_num,
                                     DTK_BLOCKED,
                                     tgt_coord.data(),
                                     tgt_num,
                                     DTK_BLOCKED,
                                     space_dim,
                                     options.c_str() );

  DTK_Map_apply( dtk_map,
                 src_field.data(),
                 DTK_INTERLEAVED,
                 tgt_field.data(),
                 DTK_INTERLEAVED,
                 field_dim );

  DTK_Map_delete( dtk_map );

}

} // end namespace
