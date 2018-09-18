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
  Init_test
  SOURCES tstInit.cpp
  COMM serial mpi
  )

TRIBITS_ADD_TEST(
  Init_test
  NAME "Init_1_test"
  ARGS "-t 1"
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )

TRIBITS_ADD_TEST(
  Init_test
  NAME "Init_2_test"
  ARGS "-t 2"
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )

TRIBITS_ADD_TEST(
  Init_test
  NAME "Init_3_test"
  ARGS "-t 3"
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )

TRIBITS_ADD_EXECUTABLE(
  C_API_test
  SOURCES tstC_API.c tstC_API.cpp
  COMM serial mpi
  ADDED_EXE_TARGET_NAME_OUT C_API_test_EXE_TARGET_NAME
  )

IF (C_API_test_EXE_TARGET_NAME)
  SET_TARGET_PROPERTIES(
    ${C_API_test_EXE_TARGET_NAME}
    PROPERTIES C_STANDARD 99
  )
ENDIF()

IF (Trilinos_ENABLE_Fortran)
    TRIBITS_ADD_EXECUTABLE(
      Fortran_API_test
      SOURCES tstFortran_API.f90 tstC_API.cpp
      COMM serial mpi
    )
ENDIF()

IF (Kokkos_ENABLE_Serial)
  TRIBITS_ADD_TEST(
    C_API_test
    NAME "C_API_Serial_test"
    ARGS "-s serial"
    STANDARD_PASS_OUTPUT
    FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
    )
  IF (Trilinos_ENABLE_Fortran)
    TRIBITS_ADD_TEST(
      Fortran_API_test
      NAME "Fortran_API_Serial_test"
      ARGS "-s serial"
      STANDARD_PASS_OUTPUT
      FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
      )
  ENDIF()
ENDIF()

IF (Kokkos_ENABLE_OpenMP)
  TRIBITS_ADD_TEST(
    C_API_test
    NAME "C_API_OpenMP_test"
    ARGS "-s openmp"
    STANDARD_PASS_OUTPUT
    FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
    )
  IF (Trilinos_ENABLE_Fortran)
    TRIBITS_ADD_TEST(
      Fortran_API_test
      NAME "Fortran_API_OpenMP_test"
      ARGS "-s openmp"
      STANDARD_PASS_OUTPUT
      FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
      )
  ENDIF()
ENDIF()

IF (Kokkos_ENABLE_Cuda)
  TRIBITS_ADD_TEST(
    C_API_test
    NAME "C_API_CUDA_test"
    ARGS "-s cuda"
    STANDARD_PASS_OUTPUT
    FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
    )
  IF (Trilinos_ENABLE_Fortran)
    TRIBITS_ADD_TEST(
      Fortran_API_test
      NAME "Fortran_API_Cuda_test"
      ARGS "-s cuda"
      STANDARD_PASS_OUTPUT
      FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
      )
  ENDIF()
ENDIF()
