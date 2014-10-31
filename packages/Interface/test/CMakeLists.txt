INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddAdvancedTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  EntityIterator_test
  SOURCES tstEntityIterator.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
