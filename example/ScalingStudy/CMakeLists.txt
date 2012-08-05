INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScalingStudy1
  SOURCES scaling_study_1.cpp
  COMM serial mpi
  DEPLIBS teuchos datatransferkit
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScalingStudy2
  SOURCES scaling_study_2.cpp
  COMM serial mpi
  DEPLIBS teuchos datatransferkit
  )
