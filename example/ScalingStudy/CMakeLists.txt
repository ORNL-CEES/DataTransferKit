INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScalingStudy1
  SOURCES scaling_study_1.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScalingStudy2
  SOURCES scaling_study_2.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScalingStudy3
  SOURCES scaling_study_3.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScalingStudy4
  SOURCES scaling_study_4.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScalingStudy5
  SOURCES scaling_study_5.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScalingStudy6
  SOURCES scaling_study_6.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )
