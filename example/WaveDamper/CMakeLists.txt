INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  WaveDamperExample
  SOURCES Wave Wave_Source Wave_Target Damper Damper_Source Damper_Target cxx_main.cpp
  COMM serial mpi  
  DEPLIBS teuchos coupler_mesh coupler_core
  )
