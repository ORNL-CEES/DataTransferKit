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
  DTK_BasicEntitySet.hpp
  DTK_Box.hpp
  DTK_BoxImpl.hpp
  DTK_EntitySet.hpp
  DTK_Entity.hpp
  DTK_EntityImpl.hpp
  DTK_Types.hpp
  DTK_MappingStatus.hpp
  DTK_Point.hpp
  DTK_Point_impl.hpp
  ) 

APPEND_SET(SOURCES
  DTK_BasicEntitySet.cpp
  DTK_Box.cpp
  DTK_BoxImpl.cpp
  DTK_EntitySet.cpp
  DTK_Entity.cpp
  DTK_EntityImpl.cpp
  DTK_MappingStatus.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_geometry
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/
)
