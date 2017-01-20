TRIBITS_ADD_EXECUTABLE_AND_TEST(
  KokkosView_test
  SOURCES tstKokkosView.cpp unit_test_main.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  KokkosHelpers_test
  SOURCES tstKokkosHelpers.cpp unit_test_main.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  PointCloud_test
  SOURCES tstPointCloud.cpp unit_test_main.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  InterpolationOperator_test
  SOURCES tstInterpolationOperator.cpp unit_test_main.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  FAIL_REGULAR_EXPRESSION "data race;leak;runtime error"
  )
