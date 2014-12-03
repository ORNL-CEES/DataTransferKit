INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  StrongScalingStudy
  SOURCES strong_scaling_study.cpp
  COMM mpi
  )

TRIBITS_ADD_EXECUTABLE(
  WeakScalingStudy
  SOURCES weak_scaling_study.cpp
  COMM mpi
  )
