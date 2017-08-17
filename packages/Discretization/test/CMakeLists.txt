# ##---------------------------------------------------------------------------##
# ## TESTS
# ##---------------------------------------------------------------------------##

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PointInCell
  SOURCES tstPointInCell.cpp unit_test_main.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )
