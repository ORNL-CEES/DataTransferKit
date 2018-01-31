# ##---------------------------------------------------------------------------##
# ## TESTS
# ##---------------------------------------------------------------------------##

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

IF (HAVE_DTK_BOOST)
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

IF (HAVE_DTK_NETCDF)
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