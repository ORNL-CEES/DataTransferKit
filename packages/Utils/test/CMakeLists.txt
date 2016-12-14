
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DBC_test
  SOURCES tstDBC.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;data race;leak;runtime error"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PredicateComposition_test
  SOURCES tstPredicateComposition.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;data race;leak;runtime error"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MpiTagConsistency_test
  SOURCES tstMpiTagConsistency.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM mpi
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;data race;leak;runtime error"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  StaticSearchTree_test
  SOURCES tstStaticSearchTree.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;data race;leak;runtime error"
  )
