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
  DTK_CommTools.hpp
  DTK_Exception.hpp
  )

# Set sources

APPEND_SET(SOURCES
  DTK_CommTools.cpp
  DTK_Exception.cpp
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
   Utils
   HEADERS ${HEADERS}
   SOURCES ${SOURCES}
   )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/)
