INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  DTK_BoundingBox.hpp
  DTK_Box.hpp
  DTK_CellTopologyFactory.hpp
  DTK_CloudSearch.hpp
  DTK_CloudSearch_impl.hpp
  DTK_CommIndexer.hpp
  DTK_CommTools.hpp
  DTK_Cylinder.hpp
  DTK_DBC.hpp
  DTK_ElementMeasure.hpp
  DTK_FieldContainer.hpp
  DTK_FieldEvaluator.hpp
  DTK_FieldIntegrator.hpp
  DTK_FieldManager.hpp
  DTK_FieldManager_def.hpp
  DTK_FieldTools.hpp
  DTK_FieldTools_def.hpp
  DTK_FieldTraits.hpp
  DTK_GeometryManager.hpp
  DTK_GeometryManager_def.hpp
  DTK_GeometryRCB.hpp
  DTK_GeometryRCB_def.hpp
  DTK_GeometryRendezvous.hpp
  DTK_GeometryRendezvous_def.hpp
  DTK_GeometryTraits.hpp
  DTK_IntegralAssemblyMap.hpp
  DTK_IntegralAssemblyMap_def.hpp
  DTK_KDTree.hpp
  DTK_KDTree_def.hpp
  DTK_MeshContainer.hpp
  DTK_MeshManager.hpp
  DTK_MeshManager_def.hpp
  DTK_MeshTools.hpp
  DTK_MeshTools_def.hpp
  DTK_MeshTraits.hpp
  DTK_MeshTraitsFieldAdapter.hpp
  DTK_MeshTypes.hpp
  DTK_Partitioner.hpp
  DTK_PartitionerFactory.hpp
  DTK_RCB.hpp
  DTK_RCB_def.hpp
  DTK_Rendezvous.hpp
  DTK_Rendezvous_def.hpp
  DTK_RendezvousMesh.hpp
  DTK_RendezvousMesh_def.hpp
  DTK_SerialPartitioner.hpp
  DTK_SharedDomainMap.hpp
  DTK_SharedDomainMap_def.hpp
  DTK_TopologyTools.hpp
  DTK_TopologyTools_def.hpp
  DTK_VolumeSourceMap.hpp
  DTK_VolumeSourceMap_def.hpp
  ) 

APPEND_SET(SOURCES
  DTK_BoundingBox.cpp
  DTK_Box.cpp
  DTK_CellTopologyFactory.cpp
  DTK_CommIndexer.cpp
  DTK_CommTools.cpp
  DTK_Cylinder.cpp
  DTK_DBC.cpp
  DTK_MeshContainer.cpp
  DTK_SerialPartitioner.cpp
  DTK_TopologyTools.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  datatransferkit
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/
)
