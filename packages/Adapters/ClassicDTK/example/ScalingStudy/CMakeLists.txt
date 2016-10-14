
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScalingStudy1
  SOURCES scaling_study_1.cpp
  COMM mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScalingStudy2
  SOURCES scaling_study_2.cpp
  COMM mpi
  )
