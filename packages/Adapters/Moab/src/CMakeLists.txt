SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  DTK_MoabEntity.hpp
  DTK_MoabEntityExtraData.hpp
  DTK_MoabEntityImpl.hpp
  DTK_MoabEntityIterator.hpp
  DTK_MoabEntityIteratorRange.hpp
  DTK_MoabEntityPredicates.hpp
  DTK_MoabEntitySet.hpp
  DTK_MoabHelpers.hpp
  DTK_MoabMeshSetIndexer.hpp
  ) 

APPEND_SET(SOURCES
  DTK_MoabEntity.cpp
  DTK_MoabEntityImpl.cpp
  DTK_MoabEntityIterator.cpp
  DTK_MoabEntityPredicates.cpp
  DTK_MoabEntitySet.cpp
  DTK_MoabHelpers.cpp
  DTK_MoabMeshSetIndexer.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_moabadapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
