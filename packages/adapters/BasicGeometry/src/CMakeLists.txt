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
  DTK_BasicGeometryEntity.hpp
  DTK_BasicGeometryEntityImpl.hpp
  DTK_BasicGeometryLocalMap.hpp
  DTK_Box.hpp
  DTK_BoxImpl.hpp
  DTK_Point.hpp
  DTK_PointImpl.hpp
  ) 

APPEND_SET(SOURCES
  DTK_BasicEntitySet.cpp
  DTK_BasicGeometryEntity.cpp
  DTK_BasicGeometryEntityImpl.cpp
  DTK_BasicGeometryLocalMap.cpp
  DTK_Box.cpp
  DTK_BoxImpl.cpp
  DTK_Point.cpp
  DTK_PointImpl.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_basicgeometryadapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/
)
