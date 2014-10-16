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
  DTK_FunctionSpace.hpp
  DTK_MapOperator.hpp
  DTK_ParallelSearch.hpp
  ) 

APPEND_SET(SOURCES
  DTK_EntityPredicates.cpp
  DTK_FunctionSpace.cpp
  DTK_MapOperator.cpp
  DTK_ParallelSearch.cpp
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
