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
