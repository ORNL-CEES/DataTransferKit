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
  DTK_EntityPredicates.hpp
  DTK_CoarseGlobalSearch.hpp
  DTK_CoarseLocalSearch.hpp
  DTK_FineLocalSearch.hpp
  ) 

APPEND_SET(SOURCES
  DTK_EntityPredicates.cpp
  DTK_CoarseGlobalSearch.cpp
  DTK_CoarseLocalSearch.cpp
  DTK_FineLocalSearch.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_operators
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/
)
