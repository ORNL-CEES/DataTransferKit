INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  StrongScalingStudy
  SOURCES strong_scaling_study.cpp
  COMM mpi
  DEPLIBS dtk_pointcloud
  )

TRIBITS_ADD_EXECUTABLE(
  WeakScalingStudy
  SOURCES weak_scaling_study.cpp
  COMM mpi
  DEPLIBS dtk_pointcloud
  )
