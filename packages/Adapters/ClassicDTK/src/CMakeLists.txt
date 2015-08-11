SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  DTK_BoundingBox.hpp
  DTK_Box.hpp
  DTK_CommIndexer.hpp
  DTK_CommTools.hpp
  DTK_Cylinder.hpp
  DTK_FieldContainer.hpp
  DTK_FieldEvaluator.hpp
  DTK_FieldManager.hpp
  DTK_FieldManager_def.hpp
  DTK_FieldTools.hpp
  DTK_FieldTools_def.hpp
  DTK_FieldTraits.hpp
  DTK_GeometryManager.hpp
  DTK_GeometryManager_def.hpp
  DTK_GeometryTraits.hpp
  DTK_MeshContainer.hpp
  DTK_MeshManager.hpp
  DTK_MeshManager_def.hpp
  DTK_MeshTools.hpp
  DTK_MeshTools_def.hpp
  DTK_MeshTraits.hpp
  DTK_MeshTraitsFieldAdapter.hpp
  DTK_MeshTypes.hpp
  DTK_SharedDomainMap.hpp
  DTK_SharedDomainMap_def.hpp
  DTK_VolumeSourceMap.hpp
  DTK_VolumeSourceMap_def.hpp
  ) 

APPEND_SET(SOURCES
  DTK_BoundingBox.cpp
  DTK_Box.cpp
  DTK_CommIndexer.cpp
  DTK_CommTools.cpp
  DTK_Cylinder.cpp
  DTK_MeshContainer.cpp
  )

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/impl)
APPEND_SET(HEADERS
  ${DIR}/DTK_ClassicGeometricEntity.hpp
  ${DIR}/DTK_ClassicGeometricEntity_impl.hpp
  ${DIR}/DTK_ClassicGeometricEntityExtraData.hpp
  ${DIR}/DTK_ClassicGeometricEntityImpl.hpp
  ${DIR}/DTK_ClassicGeometricEntityImpl_impl.hpp
  ${DIR}/DTK_ClassicGeometricEntityLocalMap.hpp
  ${DIR}/DTK_ClassicGeometricEntityLocalMap_impl.hpp
  ${DIR}/DTK_ClassicMesh.hpp
  ${DIR}/DTK_ClassicMesh_impl.hpp
  ${DIR}/DTK_ClassicMeshElement.hpp
  ${DIR}/DTK_ClassicMeshElement_impl.hpp
  ${DIR}/DTK_ClassicMeshElementImpl.hpp
  ${DIR}/DTK_ClassicMeshElementImpl_impl.hpp
  ${DIR}/DTK_ClassicMeshElementLocalMap.hpp
  ${DIR}/DTK_ClassicMeshElementLocalMap_impl.hpp
  ${DIR}/DTK_ClassicMeshEntitySet.hpp
  ${DIR}/DTK_ClassicMeshEntitySet_impl.hpp
  ${DIR}/DTK_ClassicMeshNodalShapeFunction.hpp
  ${DIR}/DTK_ClassicMeshNodalShapeFunction_impl.hpp
  ) 

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  classic_dtk
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
