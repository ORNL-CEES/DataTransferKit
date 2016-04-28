# Search tests.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CoarseGlobalSearch_test
  SOURCES tstCoarseGlobalSearch.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CoarseLocalSearch_test
  SOURCES tstCoarseLocalSearch.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  FineLocalSearch_test
  SOURCES tstFineLocalSearch.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ParallelSearch_test
  SOURCES tstParallelSearch.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# SharedDomain tests.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrixFill_test
  SOURCES tstCrsMatrixFill.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConsistentInterpolationOperator_test
  SOURCES tstConsistentInterpolationOperator.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

# PointCloud tests.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RadialBasis_test
  SOURCES tstRadialBasis.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CloudDomain_test
  SOURCES tstCloudDomain.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SplineInterpolationPairing_test
  SOURCES tstSplineInterpolationPairing.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CenterDistributor_test
  SOURCES tstCenterDistributor.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PolynomialMatrix_test
  SOURCES tstPolynomialMatrix.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SplineProlongationOperator_test
  SOURCES tstSplineProlongationOperator.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PointCloudOperators_test
  SOURCES tstPointCloudOperators.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NodeToNodeOperator_test
  SOURCES tstNodeToNodeOperator.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  PointCloudOperatorsXML
  SOURCE_FILES spline_interpolation_test_radius.xml spline_interpolation_test_knn.xml mls_test_radius.xml mls_test_knn.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS PointCloudOperators_test
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  NodeToNodeOperatorXML
  SOURCE_FILES node_to_node_test.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS NodeToNodeOperator_test
)
