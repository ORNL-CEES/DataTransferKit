INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE(
  STKConsistentInterpolation
  SOURCES consistent_interpolation.cpp
  COMM serial mpi
  )
