INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScalingStudy1
  SOURCES scaling_study_1.cpp
  COMM serial mpi
  DEPLIBS teuchos dtk_utils dtk_mesh dtk_field
  )
