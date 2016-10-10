
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  WaveDamperExample
  SOURCES Wave WaveAdapter WaveEvaluator Damper DamperAdapter DamperEvaluator cxx_main.cpp
  NUM_MPI_PROCS 1
  )

SET(NUM_MPI_PROCS_LIST 2 3 4)
FOREACH(NUM_MPI_PROCS ${NUM_MPI_PROCS_LIST})
  TRIBITS_ADD_TEST(
    WaveDamperExample
    COMM mpi
    NUM_MPI_PROCS ${NUM_MPI_PROCS}
  )
ENDFOREACH()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  WaveDamperExample2
  SOURCES Wave WaveAdapter WaveEvaluator Damper DamperAdapter DamperEvaluator cxx_main_2.cpp
  COMM serial mpi
  )
