INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddAdvancedTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RadialBasis_test
  SOURCES tstRadialBasis.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CloudDomain_test
  SOURCES tstCloudDomain.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SplineInterpolationPairing_test
  SOURCES tstSplineInterpolationPairing.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CenterDistributor_test
  SOURCES tstCenterDistributor.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SplineInterpolator_test
  SOURCES tstSplineInterpolator.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )
