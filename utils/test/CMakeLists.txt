INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddAdvancedTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CommTools_test
  SOURCES tstCommTools.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
