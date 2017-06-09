# ##---------------------------------------------------------------------------##
# ## TESTS
# ##---------------------------------------------------------------------------##

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  KokkosView_test
  SOURCES tstKokkosView.cpp unit_test_main.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  View_test
  SOURCES tstView.cpp unit_test_main.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  ${DEPS}
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  UserApplication_test
  SOURCES tstUserApplication.cpp unit_test_main.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )

TRIBITS_ADD_EXECUTABLE(
  C_API_test
  SOURCES tstC_API.c tstC_API.cpp
  COMM serial mpi
  )

IF (Kokkos_ENABLE_Serial)
  TRIBITS_ADD_TEST(
    C_API_test
    NAME "C_API_Serial_test"
    ARGS "-s serial"
    STANDARD_PASS_OUTPUT
    FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
    )
ENDIF()

IF (Kokkos_ENABLE_OpenMP)
  TRIBITS_ADD_TEST(
    C_API_test
    NAME "C_API_OpenMP_test"
    ARGS "-s openmp"
    STANDARD_PASS_OUTPUT
    FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
    )
ENDIF()

IF (Kokkos_ENABLE_Cuda)
  TRIBITS_ADD_TEST(
    C_API_test
    NAME "C_API_CUDA_test"
    ARGS "-s cuda"
    STANDARD_PASS_OUTPUT
    FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
    )
ENDIF()
