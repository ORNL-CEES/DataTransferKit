INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  GeometryToGeometryExample
  SOURCES cxx_main.cpp
  COMM serial mpi  
  DEPLIBS datatransferkit
  )
