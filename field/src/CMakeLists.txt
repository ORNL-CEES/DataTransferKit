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
  DTK_FieldEvaluator.hpp	
  DTK_FieldManager.hpp
  DTK_FieldManager_def.hpp
  DTK_FieldTools.hpp
  DTK_FieldTools_def.hpp
  DTK_FieldTraits.hpp
  DTK_FieldVersion.hpp
  DTK_MeshAdapater.hpp
  )

# Set sources

APPEND_SET(SOURCES
  DTK_FieldVersion.cpp
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
   Field
   HEADERS ${HEADERS}
   SOURCES ${SOURCES}
   )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/)
