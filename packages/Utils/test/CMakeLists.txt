TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Miscellaneous_test
  SOURCES tstMiscellaneous.cpp unit_test_main.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;data race;leak;runtime error"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DBC_test
  SOURCES tstDBC.cpp unit_test_main.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;data race;leak;runtime error"
  )
