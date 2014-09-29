
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GeometryToMeshExample
  SOURCES cxx_main.cpp
  COMM serial mpi  
  DEPLIBS datatransferkit
  )
