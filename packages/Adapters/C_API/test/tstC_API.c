#include <DTK_C_API.h>

#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  unsigned src_num = 2;
  unsigned tgt_num = 1;
  int space_dim = 2;
  int field_dim = 1;

  double* src_coord = (double*) malloc(src_num*space_dim*sizeof(double));
  double* src_field = (double*) malloc(src_num*field_dim*sizeof(double));
  double* tgt_coord = (double*) malloc(tgt_num*space_dim*sizeof(double));
  double* tgt_field = (double*) malloc(tgt_num*field_dim*sizeof(double));

  MPI_Comm comm = MPI_COMM_WORLD;
  int comm_rank;
  int comm_size;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  src_coord[0] = 0;
  src_coord[1] = comm_rank;
  src_field[0] = src_coord[1];
  src_coord[2] = 0;
  src_coord[3] = comm_rank + 0.5;
  src_field[1] = src_coord[3];
  tgt_coord[0] = 0;
  tgt_coord[1] = comm_size - comm_rank - 1;
  tgt_field[0] = 255;

  char* options = "{ \"Map Type\": \"Node To Node\" }";

  DTK_Map* dtk_map = DTK_Map_create( comm,
      src_coord, src_num, DTK_INTERLEAVED,
      tgt_coord, tgt_num, DTK_BLOCKED,
      space_dim, options );
      
  DTK_Map_apply( dtk_map,
      src_field, DTK_INTERLEAVED,
      tgt_field, DTK_INTERLEAVED,
      field_dim, false );

  DTK_Map_delete( dtk_map );

  assert( tgt_field[0] == tgt_coord[1] );

  free(src_coord);
  free(src_field);
  free(tgt_coord);
  free(tgt_field);

  MPI_Finalize();

  if (comm_rank == 0)
    printf("End Result: TEST PASSED\n");

  return 0;
}
