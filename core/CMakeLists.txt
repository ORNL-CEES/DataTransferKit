
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsSubPackageMacros)

#
# A) Package-specific configuration options
#

TRIBITS_SUBPACKAGE(Core)

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
  Coupler_Data_Source.hpp
  Coupler_Data_Target.hpp
  Coupler_Data_Field.hpp
  Coupler_Data_Field.t.hpp
  )

# Set sources

APPEND_SET(SOURCES
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
  coupler_core
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/coupler/core)

#
# D) Define the unit tests
#
TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_SUBPACKAGE_POSTPROCESS()
