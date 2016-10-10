
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# libmesh only example
TRIBITS_ADD_EXECUTABLE(
  LibmeshInterpolation
  SOURCES interpolation.cpp
  COMM serial mpi
  )

# moab/libmesh examples
IF ( TPL_ENABLE_MOAB )

  TRIBITS_ADD_EXECUTABLE(
  MoabToLibmesh2D
  SOURCES tstMoabtoLibmesh2D.cpp
  COMM mpi
  )

  TRIBITS_ADD_EXECUTABLE(
  LibmeshToMoabBlockVar2D
  SOURCES tstLibmeshtoMoabBlockVar2D.cpp
  COMM mpi
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

ENDIF()

# copy dtk data files
TRIBITS_COPY_FILES_TO_BINARY_DIR(
  LibmeshInterpolationINP
  SOURCE_FILES input.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS LibmeshInterpolation
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  LibmeshInterpolationEXO
  SOURCE_FILES cube_mesh.exo pincell_mesh.exo sahex.e sahex1_unic.h5m proteus_fumex_2d.h5m proteus_fumex_2d.exo
  SOURCE_DIR ${DTK_DATA_DIR}/exodus/
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS LibmeshInterpolation
)
