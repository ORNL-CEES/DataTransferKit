
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
  DTK_CoarseGlobalSearch.hpp
  DTK_CoarseLocalSearch.hpp
  DTK_ConsistentInterpolationOperator.hpp
  DTK_ConsistentInterpolationOperator_impl.hpp
  DTK_FineLocalSearch.hpp
  DTK_ParallelSearch.hpp
  ) 

APPEND_SET(SOURCES
  DTK_CoarseGlobalSearch.cpp
  DTK_CoarseLocalSearch.cpp
  DTK_FineLocalSearch.cpp
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
