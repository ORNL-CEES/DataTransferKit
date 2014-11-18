INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  STKConsistentInterpolation
  SOURCES consistent_interpolation.cpp
  COMM serial mpi
  )
