INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddAdvancedTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CoarseLocalSearch_test
  SOURCES tstCoarseLocalSearch.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ParallelSearch_test
  SOURCES tstParallelSearch.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
