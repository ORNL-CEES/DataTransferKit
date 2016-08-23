##---------------------------------------------------------------------------##
# Hex mesh reference library for unit testing
##---------------------------------------------------------------------------##
SET(SOURCES
  reference_implementation/DTK_ReferenceHex.cpp
  reference_implementation/DTK_ReferenceHexImpl.cpp
  reference_implementation/DTK_ReferenceHexIntegrationRule.cpp
  reference_implementation/DTK_ReferenceHexLocalMap.cpp
  reference_implementation/DTK_ReferenceHexMesh.cpp
  reference_implementation/DTK_ReferenceHexShapeFunction.cpp
  reference_implementation/DTK_ReferenceNode.cpp
  reference_implementation/DTK_ReferenceNodeImpl.cpp
  )

TRIBITS_ADD_LIBRARY(
  dtk_hex_test_reference
  SOURCES ${SOURCES}
  TESTONLY
  )

##---------------------------------------------------------------------------##
# Reference hex implementation tests
##---------------------------------------------------------------------------##
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ReferenceNode_test
  SOURCES tstReferenceNode.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  TESTONLYLIBS dtk_hex_test_reference
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ReferenceHex_test
  SOURCES tstReferenceHex.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  TESTONLYLIBS dtk_hex_test_reference  
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ReferenceHexLocalMap_test
  SOURCES tstReferenceHexLocalMap.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  TESTONLYLIBS dtk_hex_test_reference  
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ReferenceHexShapeFunction_test
  SOURCES tstReferenceHexShapeFunction.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  TESTONLYLIBS dtk_hex_test_reference  
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ReferenceHexIntegrationRule_test
  SOURCES tstReferenceHexIntegrationRule.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  TESTONLYLIBS dtk_hex_test_reference  
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ReferenceHexMesh_test
  SOURCES tstReferenceHexMesh.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  TESTONLYLIBS dtk_hex_test_reference  
  )

##---------------------------------------------------------------------------##
# Search tests.
##---------------------------------------------------------------------------##
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

##---------------------------------------------------------------------------##
# SharedDomain tests.
##---------------------------------------------------------------------------##
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

##---------------------------------------------------------------------------##
# PointCloud tests.
##---------------------------------------------------------------------------##
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
  NodeToNodeOperator_test
  SOURCES tstNodeToNodeOperator.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
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
  VirtualWork_test
  SOURCES tstVirtualWork.cpp ${TEUCHOS_STD_PARALLEL_UNIT_TEST_MAIN}
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  PointCloudOperatorsXML
  SOURCE_FILES spline_interpolation_test_radius.xml spline_interpolation_test_knn.xml mls_test_radius.xml mls_test_knn.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS PointCloudOperators_test VirtualWork_test
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  NodeToNodeOperatorXML
  SOURCE_FILES node_to_node_test.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS NodeToNodeOperator_test
)
