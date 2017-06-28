
#
# A) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
IF (Trilinos_ENABLE_Fortran)
    APPEND_GLOB(SOURCES ${DIR}/DTK_Fortran_API.F90)
    APPEND_GLOB(SOURCES ${DIR}/DTK_Fortran_wrap.cpp)
ENDIF()

TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# B) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_interface
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  ADDED_LIB_TARGET_NAME_OUT DTK_INTERFACE_LIBNAME
  )

IF (Trilinos_ENABLE_Fortran)
  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/datatransferkit.mod DESTINATION include)
ENDIF()

# We need to set the linker language explicitly here for CUDA builds.
SET_PROPERTY(
  TARGET ${DTK_INTERFACE_LIBNAME}
  APPEND PROPERTY LINKER_LANGUAGE CXX
  )
