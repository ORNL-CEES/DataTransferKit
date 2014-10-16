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
  DTK_Entity.hpp
  DTK_EntityExtraData.hpp
  DTK_EntityImpl.hpp
  DTK_EntityIterator.hpp
  DTK_EntityLocalMap.hpp
  DTK_EntitySet.hpp
  DTK_EntityShapeFunction.hpp
  DTK_Types.hpp
  ) 

APPEND_SET(SOURCES
  DTK_Entity.cpp
  DTK_EntityImpl.cpp
  DTK_EntityIterator.cpp
  DTK_EntityLocalMap.cpp
  DTK_EntitySet.cpp
  DTK_EntityShapeFunction.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_interface
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/
)
