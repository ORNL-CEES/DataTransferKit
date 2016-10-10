SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  DTK_LibmeshAdjacencies.hpp
  DTK_LibmeshEntity.hpp
  DTK_LibmeshEntityExtraData.hpp
  DTK_LibmeshEntityImpl.hpp
  DTK_LibmeshEntityIntegrationRule.hpp
  DTK_LibmeshEntityIterator.hpp
  DTK_LibmeshEntityIterator_impl.hpp
  DTK_LibmeshEntityLocalMap.hpp
  DTK_LibmeshEntitySet.hpp
  DTK_LibmeshHelpers.hpp
  DTK_LibmeshManager.hpp
  DTK_LibmeshNodalShapeFunction.hpp
  DTK_LibmeshVariableField.hpp
  )

APPEND_SET(SOURCES
  DTK_LibmeshAdjacencies.cpp
  DTK_LibmeshEntity.cpp
  DTK_LibmeshEntityImpl.cpp
  DTK_LibmeshEntityIntegrationRule.cpp
  DTK_LibmeshEntityLocalMap.cpp
  DTK_LibmeshEntitySet.cpp
  DTK_LibmeshManager.cpp
  DTK_LibmeshNodalShapeFunction.cpp
  DTK_LibmeshVariableField.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_libmeshadapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
