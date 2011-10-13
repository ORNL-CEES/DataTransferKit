# CMakeLists.txt for dtransfer/src/core
# INCLUDE(PackageLibraryMacros)
# INCLUDE(PackageMacros)

# Package-specific configuration
# SET(DTRANSFER_DEBUG ${HAVE_DTRANSFER_DEBUG})
# PACKAGE_CONFIGURE_FILE(config.h)

# Tests
#PACKAGE_ADD_TEST_DIRECTORIES(test)

# Headers and sources
SET(HEADERS
  Data_Transfer_Manager.hh
  Data_Transfer_Manager.i.hh
  Data_Transfer_Manager.t.hh
  Transfer_Evaluator.hh
  Transfer_Evaluator.i.hh
  Transfer_Evaluator.t.hh
  Transfer_Map.hh
  Transfer_Map.i.hh
  Transfer_Map.t.hh
)

SET(SOURCES
  Data_Transfer_Manager.cc
  Transfer_Evaluator.cc
  Transfer_Map.cc
)

# Targets and install

# PACKAGE_ADD_LIBRARY(
#   ${DTRANSFER_PREFIX}core
#   SOURCES ${SOURCES}
# )

INSTALL(FILES
  ${HEADERS}
  DESTINATION
  ${CMAKE_INSTALL_PREFIX}/include/core
  )

# Library
ADD_LIBRARY( core ${SOURCES} )
INSTALL_TARGETS( ${CMAKE_INSTALL_PREFIX}/lib core )
