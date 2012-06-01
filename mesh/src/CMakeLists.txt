# CMakeLists.txt for src

INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set headers

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  DTK_CellTopologyFactory.hpp
  DTK_KDTree.hpp
  DTK_KDTree_def.hpp
  DTK_Mesh.hpp
  DTK_Mesh_def.hpp
  DTK_TopologyTools.hpp
  )

# Set sources

APPEND_SET(SOURCES
  DTK_CellTopologyFactory.cpp
  DTK_TopologyTools.cpp
  )

# Set addition headers and source based on package-specific
# configuration options

IF(${PACKAGE_NAME}_ENABLE_EXTENDED)
APPEND_SET(HEADERS
    )

  APPEND_SET(SOURCES
    )

ENDIF()

#
# C) Define the targets for install
#

TRIBITS_ADD_LIBRARY(
   Mesh
   HEADERS ${HEADERS}
   SOURCES ${SOURCES}
   )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/)
