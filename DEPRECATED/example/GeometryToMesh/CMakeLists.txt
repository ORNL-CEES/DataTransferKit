INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  GeometryToMeshExample
  SOURCES cxx_main.cpp
  COMM serial mpi  
  DEPLIBS datatransferkit
  )
