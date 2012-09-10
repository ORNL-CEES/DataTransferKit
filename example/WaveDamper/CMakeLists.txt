INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  WaveDamperExample
  SOURCES Wave WaveAdapter WaveEvaluator Damper DamperAdapter DamperEvaluator cxx_main.cpp
  COMM mpi  
  DEPLIBS teuchos datatransferkit
  )
