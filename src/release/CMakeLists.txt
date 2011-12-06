
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsSubPackageMacros)

#
# A) Package-specific configuration options
#

TRIBITS_SUBPACKAGE(Release)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# NOTE: Above, we add the XXX_config.h file so tht MS Visual Studio
# will know to include this file in the dependency analysis

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  config.h.in
  Release.hh
  )

# Set sources

APPEND_SET(SOURCES
  Release.cc
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
  release
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DCOUPLER_LIB_EXPORTS_MODE
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/release)

TRIBITS_SUBPACKAGE_POSTPROCESS()