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
  Point.hh
  Point.i.hh
  Point.t.hh
)

SET(SOURCES
  Point.cc
)

# Targets and install

# PACKAGE_ADD_LIBRARY(
#   ${DTRANSFER_PREFIX}core
#   SOURCES ${SOURCES}
# )

INSTALL(FILES
  ${HEADERS}
  DESTINATION
  ${CMAKE_INSTALL_PREFIX}/include/mesh
  )

# Library
ADD_LIBRARY( mesh ${SOURCES} )
INSTALL_TARGETS( ${CMAKE_INSTALL_PREFIX}/lib mesh )
