INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  WaveDamperExample
  SOURCES Wave WaveAdapter WaveEvaluator Damper DamperAdapter DamperEvaluator cxx_main.cpp
  COMM serial mpi  
  DEPLIBS datatransferkit
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  WaveDamperExample2
  SOURCES Wave WaveAdapter WaveEvaluator Damper DamperAdapter DamperEvaluator cxx_main_2.cpp
  COMM serial mpi  
  DEPLIBS datatransferkit
  )
