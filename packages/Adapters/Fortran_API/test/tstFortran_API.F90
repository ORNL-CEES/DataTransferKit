program dtk_fortran_api

  use DataTransferKit
  use iso_c_binding
  implicit none

  include 'mpif.h'

  integer :: ierr, comm_rank, comm_size
  logical :: test
  type(c_ptr) :: dtk_map
  real(c_double), allocatable, target :: src_coord(:), src_field(:)
  integer(c_size_t) :: src_num
  real(c_double), allocatable, target :: tgt_coord(:), tgt_field(:)
  integer(c_size_t) :: tgt_num
  integer(c_int) :: space_dim, field_dim
  character(len=256) :: options
  logical(c_bool) :: apply_transpose

  call MPI_INIT(ierr)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, comm_rank, ierr)

  src_num = 2
  tgt_num = 1
  space_dim = 2
  field_dim = 1
  allocate(src_coord(src_num*space_dim))
  allocate(src_field(src_num*field_dim))
  allocate(tgt_coord(tgt_num*space_dim))
  allocate(tgt_field(tgt_num*field_dim))
  src_coord(1) = 0
  src_coord(2) = comm_rank
  src_field(1) = src_coord(2)
  src_coord(3) = 0
  src_coord(4) = comm_rank + 0.5
  src_field(2) = src_coord(4)
  tgt_coord(1) = 0
  tgt_coord(2) = comm_size - comm_rank - 1
  tgt_field(1) = 255
  test = .FALSE.

  apply_transpose = .FALSE.
  options = '{ "Map Type": "Node To Node" }'//c_null_char
  
  dtk_map = DTK_Map_create(MPI_COMM_WORLD, &
      c_loc(src_coord), src_num, DTK_INTERLEAVED, &
      c_loc(tgt_coord), tgt_num, DTK_BLOCKED, &
      space_dim, options)

  call DTK_Map_apply(dtk_map, c_loc(src_field), DTK_INTERLEAVED, &
      c_loc(tgt_field), DTK_INTERLEAVED, field_dim, apply_transpose)

  call MPI_ALLREDUCE(tgt_field(1) == tgt_coord(2), test, 1, MPI_LOGICAL, &
      MPI_LAND, MPI_COMM_WORLD, ierr)
  if (comm_rank == 0) then
    print *, "Does this test pass [T/F]:", test
  end if

  call DTK_Map_delete(dtk_map)

  call MPI_FINALIZE(ierr)

  if (comm_rank == 0) then
    print *, "End Result: TEST PASSED"
  end if

end program
