##---------------------------------------------------------------------------##
## coupler/CMakeLists.txt
## Thomas M. Evans
## Wednesday June 22 11:14:40 2011
##---------------------------------------------------------------------------##
## Copyright (C) 2011 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## CMAKE for coupler
##---------------------------------------------------------------------------##

INCLUDE(PackageLibraryMacros)
INCLUDE(NemesisPackageMacros)
INCLUDE(PackageMacros)

##---------------------------------------------------------------------------##
## PACKAGE-SPECIFIC OPTIONS
##---------------------------------------------------------------------------##


##---------------------------------------------------------------------------##
## PACKAGE CONFIGURE FILE
##---------------------------------------------------------------------------##

NEMESIS_PACKAGE_CONFIGURE_FILE(config.h)

##---------------------------------------------------------------------------##
## INCLUDES
##---------------------------------------------------------------------------##

INCLUDE_DIRECTORIES(${NEMESIS_BINARY_SRC_PATH}
  ${NEMESIS_SRC_PATH})

##---------------------------------------------------------------------------##
## HEADER AND SOURCE FILES
##---------------------------------------------------------------------------##

# initial declaration
SET(HEADERS "")
SET(NOINSTALLHEADERS "")
SET(SOURCES "")

# C/C++ headers
FILE(GLOB HEADERS *.hh)

# C/C++ source
FILE(GLOB SOURCES *.cc)

# configure file
APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/config.h)

##---------------------------------------------------------------------------##
## PACKAGE TARGETS AND INSTALL
##---------------------------------------------------------------------------##

PACKAGE_ADD_LIBRARY(
  ${DENOVO_PREFIX}coupler
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES})

INSTALL(FILES
  ${HEADERS}
  DESTINATION
  ${CMAKE_INSTALL_PREFIX}/include/coupler)

##---------------------------------------------------------------------------##
##                   end of coupler/CMakeLists.txt
##---------------------------------------------------------------------------##
