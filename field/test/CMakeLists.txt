INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddAdvancedTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConsistentInterpolation_test1
  SOURCES tstConsistentInterpolation1.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConsistentInterpolation_test2
  SOURCES tstConsistentInterpolation2.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConsistentInterpolation_test3
  SOURCES tstConsistentInterpolation3.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConsistentInterpolation_test4
  SOURCES tstConsistentInterpolation4.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TransferOperator_test
  SOURCES tstTransferOperator.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  )
