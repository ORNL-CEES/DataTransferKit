
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsSubPackageMacros)

#
# A) Package-specific configuration options
#

TRIBITS_SUBPACKAGE(core)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

# NOTE: Above, we add the XXX_config.h file so tht MS Visual Studio
# will know to include this file in the dependency analysis

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Packing_Utils.hh
  Packing_Utils.i.hh
  Message_Buffer.hh
  LG_Indexer.hh
  LG_Indexer.i.hh
  Transfer_Map.hh
  Transfer_Map.i.hh
  Transfer_Data_Source.hh
  Transfer_Data_Target.hh
  Transfer_Data_Field.hh
  Mapper.hh
  Messenger.hh
  Data_Transfer_Manager.hh
  )

# Set sources

APPEND_SET(SOURCES
  Transfer_Data_Field.cc
  Mapper.cc
  Messenger.cc
  Data_Transfer_Manager.cc
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
  core
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DCOUPLER_LIB_EXPORTS_MODE
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/core)

#
# D) Define the unit tests
#
TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()