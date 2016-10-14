
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  STKParallelSearch
  SOURCES parallel_search.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  STKInterpolation
  SOURCES interpolation.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  STKInlineInterpolation
  SOURCES inline_interpolation.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  STKRepeatedTransfer
  SOURCES repeated_transfer.cpp
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  STKInterpolationINP
  SOURCE_FILES input.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS STKInterpolation STKRepeatedTransfer STKInlineInterpolation
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  STKInterpolationEXO
  SOURCE_FILES cube_mesh.exo pincell_mesh.exo
  SOURCE_DIR ${DTK_DATA_DIR}/exodus/
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS STKInterpolation STKRepeatedTransfer STKInlineInterpolation
)
