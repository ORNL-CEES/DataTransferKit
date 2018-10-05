SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)


# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_mapfactory
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS dtk_utils dtk_interface dtk_meshfree
  ADDED_LIB_TARGET_NAME_OUT DTK_MAPFACTORY_LIBNAME
  )

# We need to set the linker language explicitly here for CUDA builds.
SET_PROPERTY(
  TARGET ${DTK_MAPFACTORY_LIBNAME}
  APPEND PROPERTY LINKER_LANGUAGE CXX
  )
