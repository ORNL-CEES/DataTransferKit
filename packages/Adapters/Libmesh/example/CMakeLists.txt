
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE(
  LibmeshInterpolation
  SOURCES interpolation.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  MoabToLibmeshParallel
  SOURCES tstMoabtoLibmeshParallelTransfer.cpp
  COMM mpi
  )

TRIBITS_ADD_EXECUTABLE(
  LibmeshToMoabParallel
  SOURCES tstLibmeshToMoabParallelTransfer.cpp
  COMM mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  LibmeshInterpolationINP
  SOURCE_FILES input.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS LibmeshInterpolation MoabToLibmeshParallel LibmeshToMoabParallel
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  LibmeshInterpolationEXO 
  SOURCE_FILES cube_mesh.exo pincell_mesh.exo sahex.e sahex1_unic.h5m
  SOURCE_DIR ${DTK_DATA_DIR}/exodus/
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS LibmeshInterpolation MoabToLibmeshParallel LibmeshToMoabParallel
)
