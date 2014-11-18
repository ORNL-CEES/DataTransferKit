INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  STKConsistentInterpolation
  SOURCES consistent_interpolation.cpp
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  STKConsistentInterpolationINP
  SOURCE_FILES input.xml ${DTKDATA_DIR}/exodus/cube_mesh.exo ${DTKDATA_DIR}/exodus/pincell_mesh.exo
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS STKConsistentInterpolation
)
