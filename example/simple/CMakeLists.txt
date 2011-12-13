INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Simple_example
  COMM serial mpi  
  NUM_MPI_PROCS 4
  SOURCES Physics_A.cpp Physics_B.cpp Source_A.cpp Target_B.cpp cxx_main.cpp
  )
