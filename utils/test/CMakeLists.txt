INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddAdvancedTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  InverseComm_test
  SOURCES tstInverseComm.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  )
