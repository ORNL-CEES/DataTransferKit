##---------------------------------------------------------------------------##
## TESTS
##---------------------------------------------------------------------------##
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NearestNeighborOperator
  SOURCES tstNearestNeighborOperator.cpp unit_test_main.cpp
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DetailsCommunicationHelpers
  SOURCES tstDetailsCommunicationHelpers.cpp unit_test_main.cpp
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )

IF (HAVE_DTK_BOOST AND ((BOOST_VERSION VERSION_EQUAL 1.62.0) OR (BOOST_VERSION VERSION_GREATER 1.62.0)))
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    CompactlySupportedRadialBasisFunctions
    SOURCES tstCompactlySupportedRadialBasisFunctions.cpp unit_test_main.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
    )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultivariatePolynomialBasis
  SOURCES tstMultivariatePolynomialBasis.cpp unit_test_main.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )

IF (HAVE_DTK_NETCDF AND DTK_DATA_DIR)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    NearestNeighborExodusGenerator
    SOURCES tstNearestNeighborExodusGenerator.cpp unit_test_main.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
    ENVIRONMENT CUDA_LAUNCH_BLOCKING=1
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    ExodusGeneratorFiles
    SOURCE_FILES coarse_sphere.exo fine_sphere.exo
    SOURCE_DIR ${DTK_DATA_DIR}/exodus/
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    EXEDEPS NearestNeighborExodusGenerator
    )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MeshfreeOperators
  SOURCES tstMeshfreeOperators.cpp unit_test_main.cpp
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MeshfreeOperatorsSimpleProblem
  SOURCES tstMeshfreeOperatorsSimpleProblem.cpp unit_test_main.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SVD
  SOURCES tstSVD.cpp unit_test_main.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )
