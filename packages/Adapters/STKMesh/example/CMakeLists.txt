INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  STKConsistentInterpolation
  SOURCES consistent_interpolation.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  STKMLSReconstruction
  SOURCES mls_reconstruction.cpp
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  STKConsistentInterpolationINP
  SOURCE_FILES input.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS STKConsistentInterpolation STKMLSReconstruction
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  STKConsistentInterpolationEXO
  SOURCE_FILES cube_mesh.exo pincell_mesh.exo
  SOURCE_DIR ${DTK_DATA_DIR}/exodus/
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS STKConsistentInterpolation STKMLSReconstruction
)
