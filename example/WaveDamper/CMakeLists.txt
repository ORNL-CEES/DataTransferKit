INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  WaveDamperExample
  SOURCES Wave WaveAdapter Damper DamperAdapter cxx_main.cpp
  COMM serial mpi  
  DEPLIBS teuchos datatransferkit
  )
