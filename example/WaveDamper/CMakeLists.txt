INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  WaveDamperExample
  SOURCES Wave Wave_Source Wave_Target Damper Damper_Source Damper_Target cxx_main.cpp
  ARGS
  COMM serial mpi  
  STANDARD_PASS_OUTPUT
  )
