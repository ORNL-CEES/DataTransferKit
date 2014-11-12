INCLUDE(TribitsLibraryMacros)

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
  DTK_IntrepidCell.hpp
  DTK_IntrepidSideCell.hpp
  DTK_STKMeshEntity.hpp
  DTK_STKMeshEntityExtraData.hpp
  DTK_STKMeshEntityImpl.hpp
  DTK_STKMeshEntityIterator.hpp
  DTK_STKMeshEntityIteratorRange.hpp
  DTK_STKMeshEntityPredicates.hpp
  DTK_STKMeshEntitySet.hpp
  ) 

APPEND_SET(SOURCES
  DTK_IntrepidCell.cpp
  DTK_IntrepidSideCell.cpp
  DTK_STKMeshEntity.cpp
  DTK_STKMeshEntityImpl.cpp
  DTK_STKMeshEntityIterator.cpp
  DTK_STKMeshEntityPredicates.cpp
  DTK_STKMeshEntitySet.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_stkmeshadapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/
)
