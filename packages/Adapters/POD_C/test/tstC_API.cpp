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
  int const comm_rank = teuchos_comm->getRank();
  int const comm_size = teuchos_comm->getSize();

  // set options using JSON format
  std::string const options = "{ "\
    "\"Map Type\": \"Moving Least Square Reconstruction\", "\
    "\"Basis Type\": \"Wendland\", "\
    "\"Basis Order\": 2, "\
    "\"RBF Radius\": 0.3 }";

  // TODO: fill these
  std::srand(123*comm_rank);
  double coord[3];
  int const space_dim = 3;
  int const field_dim = 2;
  unsigned const src_num = 3000;
  std::vector<double> src_coord(space_dim*src_num);
  std::vector<double> src_field(field_dim*src_num);
  for (unsigned i = 0; i < src_num; ++i)
  {
    coord[0] = (double) std::rand() / (double) RAND_MAX + comm_rank;
    coord[1] = (double) std::rand() / (double) RAND_MAX;
    coord[2] = (double) std::rand() / (double) RAND_MAX;
    // source coordinates data blocked
    src_coord[i+0*src_num] = coord[0];
    src_coord[i+1*src_num] = coord[1];
    src_coord[i+2*src_num] = coord[2];
    // source field data interleaved
    src_field[field_dim*i+0] = coord[0];
    src_field[field_dim*i+1] = coord[2];
  }
  unsigned const tgt_num = 50;
  std::vector<double> tgt_coord(space_dim*tgt_num);
  std::vector<double> tgt_field(field_dim*tgt_num);
  std::vector<double> gold_data(field_dim*tgt_num);
  for (unsigned i = 0; i < tgt_num; ++i)
  {
    coord[0] = (double) std::rand() / (double) RAND_MAX + comm_size - comm_rank - 1;
    coord[1] = (double) std::rand() / (double) RAND_MAX;
    coord[2] = (double) std::rand() / (double) RAND_MAX;
    // target coordinates data interleaved
    tgt_coord[space_dim*i+0] = coord[0];
    tgt_coord[space_dim*i+1] = coord[1];
    tgt_coord[space_dim*i+2] = coord[2];
    // target field data blocked
    gold_data[i+tgt_num*0] = coord[0];
    gold_data[i+tgt_num*1] = coord[2];
  }
  TEST_EQUALITY( src_coord.size(), space_dim * src_num );
  TEST_EQUALITY( tgt_coord.size(), space_dim * tgt_num );

  DTK_Map* dtk_map = DTK_Map_create( mpi_comm,
                                     src_coord.data(),
                                     src_num,
                                     DTK_BLOCKED,
                                     tgt_coord.data(),
                                     tgt_num,
                                     DTK_INTERLEAVED,
                                     space_dim,
                                     options.c_str() );

  DTK_Map_apply( dtk_map,
                 src_field.data(),
                 DTK_INTERLEAVED,
                 tgt_field.data(),
                 DTK_BLOCKED,
                 field_dim );

  DTK_Map_delete( dtk_map );

  double const rel_tol = 1e-8;
  double const abs_tol = 1e-6;
  for (unsigned i = 0; i < field_dim*tgt_num; ++i)
  {
    TEST_FLOATING_EQUALITY( gold_data[i], tgt_field[i], rel_tol );
    TEST_COMPARE( std::abs(gold_data[i] - tgt_field[i]), <, abs_tol );
  }
}

} // end namespace
