# CMakeLists.txt for src

INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

TRIBITS_ADD_OPTION_AND_DEFINE(
	DataTransferKit_ENABLE_DBC
	HAVE_DTK_DBC
	"Enable Design-by-Contract checks. Warning, use for debug checking but disable for production runs as it incurs a significant overhead."
	OFF
)

TRIBITS_CONFIGURE_FILE(DataTransferKit_config.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  DTK_CommIndexer.hpp
  DTK_CommTools.hpp
  DTK_Assertion.hpp
  )

# Set sources

APPEND_SET(SOURCES
  DTK_CommIndexer.cpp
  DTK_CommTools.cpp
  DTK_Assertion.cpp
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
