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
  DTK_Assertion.hpp
  DTK_BoundingBox.hpp
  DTK_CellTopologyFactory.hpp
  DTK_CommIndexer.hpp
  DTK_CommTools.hpp
  DTK_FieldContainer.hpp
  DTK_FieldEvaluator.hpp
  DTK_FieldManager.hpp
  DTK_FieldManager_def.hpp
  DTK_FieldTools.hpp
  DTK_FieldTools_def.hpp
  DTK_FieldTraits.hpp
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
  DTK_RCB.hpp
  DTK_RCB_def.hpp
  DTK_Rendezvous.hpp
  DTK_Rendezvous_def.hpp
  DTK_RendezvousMesh.hpp
  DTK_RendezvousMesh_def.hpp
  DTK_SharedDomainMap.hpp
  DTK_SharedDomainMap_def.hpp
  DTK_TopologyTools.hpp
  ) 

APPEND_SET(SOURCES
  DTK_Assertion.cpp
  DTK_BoundingBox.cpp
  DTK_CellTopologyFactory.cpp
  DTK_CommIndexer.cpp
  DTK_CommTools.cpp
  DTK_KDTree.cpp
  DTK_MeshContainer.cpp
  DTK_RendezvousMesh.cpp
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
