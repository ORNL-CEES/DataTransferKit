INCLUDE(TribitsAddExecutableAndTest)

TRIBITS_ADD_EXECUTABLE(
  WaveDamperExample
  SOURCES Wave WaveAdapter WaveEvaluator Damper DamperAdapter DamperEvaluator cxx_main.cpp
  )

TRIBITS_ADD_EXECUTABLE(
  WaveDamperExample2
  SOURCES Wave WaveAdapter WaveEvaluator Damper DamperAdapter DamperEvaluator cxx_main_2.cpp
  )
