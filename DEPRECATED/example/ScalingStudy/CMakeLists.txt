INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  ScalingStudy1
  SOURCES scaling_study_1.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )

TRIBITS_ADD_EXECUTABLE(
  ScalingStudy2
  SOURCES scaling_study_2.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )

TRIBITS_ADD_EXECUTABLE(
  ScalingStudy3
  SOURCES scaling_study_3.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )

TRIBITS_ADD_EXECUTABLE(
  ScalingStudy4
  SOURCES scaling_study_4.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )

TRIBITS_ADD_EXECUTABLE(
  ScalingStudy5
  SOURCES scaling_study_5.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )

TRIBITS_ADD_EXECUTABLE(
  ScalingStudy6
  SOURCES scaling_study_6.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )

TRIBITS_ADD_EXECUTABLE(
  ScalingStudy7
  SOURCES scaling_study_7.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )

TRIBITS_ADD_EXECUTABLE(
  ScalingStudy8
  SOURCES scaling_study_8.cpp
  COMM mpi
  DEPLIBS datatransferkit
  )
