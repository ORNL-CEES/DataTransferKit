INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  peaks_2d
  SOURCES peaks_2d.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  )
