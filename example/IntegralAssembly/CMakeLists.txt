INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  IntegralAssemblyExample
  SOURCES cxx_main.cpp
  COMM serial mpi  
  DEPLIBS datatransferkit
  )
