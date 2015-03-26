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
  DTK_MoabEntityLocalMap.hpp
  DTK_MoabEntityPredicates.hpp
  DTK_MoabEntitySet.hpp
  DTK_MoabHelpers.hpp
  DTK_MoabManager.hpp
  DTK_MoabMeshSetIndexer.hpp
  DTK_MoabNodalShapeFunction.hpp
  DTK_MoabTagField.hpp
  DTK_MoabTagField_impl.hpp
  ) 

APPEND_SET(SOURCES
  DTK_MoabEntity.cpp
  DTK_MoabEntityImpl.cpp
  DTK_MoabEntityIterator.cpp
  DTK_MoabEntityLocalMap.cpp
  DTK_MoabEntityPredicates.cpp
  DTK_MoabEntitySet.cpp
  DTK_MoabHelpers.cpp
  DTK_MoabManager.cpp
  DTK_MoabMeshSetIndexer.cpp
  DTK_MoabNodalShapeFunction.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_moabadapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
