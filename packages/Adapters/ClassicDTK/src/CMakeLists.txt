SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  DTK_Classic_BoundingBox.hpp
  DTK_Classic_Box.hpp
  DTK_Classic_Cylinder.hpp
  DTK_Classic_FieldContainer.hpp
  DTK_Classic_FieldEvaluator.hpp
  DTK_Classic_FieldManager.hpp
  DTK_Classic_FieldManager_def.hpp
  DTK_Classic_FieldTools.hpp
  DTK_Classic_FieldTools_def.hpp
  DTK_Classic_FieldTraits.hpp
  DTK_Classic_GeometryManager.hpp
  DTK_Classic_GeometryManager_def.hpp
  DTK_Classic_GeometryTraits.hpp
  DTK_Classic_MeshContainer.hpp
  DTK_Classic_MeshManager.hpp
  DTK_Classic_MeshManager_def.hpp
  DTK_Classic_MeshTools.hpp
  DTK_Classic_MeshTools_def.hpp
  DTK_Classic_MeshTraits.hpp
  DTK_Classic_MeshTraitsFieldAdapter.hpp
  DTK_Classic_MeshTypes.hpp
  DTK_Classic_SharedDomainMap.hpp
  DTK_Classic_SharedDomainMap_def.hpp
  DTK_Classic_VolumeSourceMap.hpp
  DTK_Classic_VolumeSourceMap_def.hpp
  ) 

APPEND_SET(SOURCES
  DTK_Classic_BoundingBox.cpp
  DTK_Classic_Box.cpp
  DTK_Classic_Cylinder.cpp
  DTK_Classic_MeshContainer.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  classic_dtk
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
