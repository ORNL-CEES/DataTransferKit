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

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  DataTransferKit_CommIndexer.hpp
  DataTransferKit_CommIndexer_Def.hpp
  DataTransferKit_CopyOperator.hpp
  DataTransferKit_CopyOperator_Def.hpp
  DataTransferKit_DataSource.hpp
  DataTransferKit_DataSource_Def.hpp
  DataTransferKit_DataTarget.hpp
  DataTransferKit_Exception.hpp
  DataTransferKit_Point.hpp
  DataTransferKit_SerializationTraits.hpp
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
  data_transfer_kit
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/)
