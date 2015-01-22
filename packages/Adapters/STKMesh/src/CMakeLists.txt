
#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  DTK_STKMeshEntity.hpp
  DTK_STKMeshEntityExtraData.hpp
  DTK_STKMeshEntityImpl.hpp
  DTK_STKMeshEntityIterator.hpp
  DTK_STKMeshEntityIteratorRange.hpp
  DTK_STKMeshEntityLocalMap.hpp
  DTK_STKMeshFieldVector.hpp
  DTK_STKMeshFieldVector_impl.hpp
  DTK_STKMeshManager.hpp
  DTK_STKMeshNodalShapeFunction.hpp
  DTK_STKMeshEntityPredicates.hpp
  DTK_STKMeshEntitySet.hpp
  DTK_STKMeshHelpers.hpp
  DTK_STKMeshHelpers_impl.hpp
  ) 

APPEND_SET(SOURCES
  DTK_STKMeshEntity.cpp
  DTK_STKMeshEntityImpl.cpp
  DTK_STKMeshEntityIterator.cpp
  DTK_STKMeshEntityLocalMap.cpp
  DTK_STKMeshManager.cpp
  DTK_STKMeshNodalShapeFunction.cpp
  DTK_STKMeshEntityPredicates.cpp
  DTK_STKMeshEntitySet.cpp
  DTK_STKMeshHelpers.cpp  
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_stkmeshadapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
