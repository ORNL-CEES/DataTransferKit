ADD_SUBDIRECTORY(Utils)

IF(DataTransferKit_ARBORX_TPL)
  # Create a dummy library
  SET(HEADERS "")
  SET(SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/dummy.cc")

  TRIBITS_ADD_LIBRARY(
    dtk_search
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    ADDED_LIB_TARGET_NAME_OUT DTK_SEARCH_LIBNAME
    )
ELSE()
  SET(ARBORX_ENABLE_MPI TRUE)
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Search/src/ArborX_Config.hpp.in
                 ${CMAKE_CURRENT_BINARY_DIR}/Search/src/ArborX_Config.hpp)

  SET(HEADERS "")
  SET(SOURCES "")

  FOREACH(dir ${CMAKE_CURRENT_SOURCE_DIR}/Search/src
              ${CMAKE_CURRENT_SOURCE_DIR}/Search/src/details
              ${CMAKE_CURRENT_SOURCE_DIR}/Search/src/geometry
              ${CMAKE_CURRENT_SOURCE_DIR}/Search/src/kokkos_ext
              ${CMAKE_CURRENT_BINARY_DIR}/Search/src )
    APPEND_GLOB(HEADERS ${dir}/*.hpp)
    APPEND_GLOB(SOURCES ${dir}/*.cpp)
    TRIBITS_INCLUDE_DIRECTORIES(${dir})
  ENDFOREACH()

  TRIBITS_ADD_LIBRARY(
    dtk_search
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS dtk_utils
    ADDED_LIB_TARGET_NAME_OUT DTK_SEARCH_LIBNAME
    )

  # We need to set the linker language explicitly here for CUDA builds.
  SET_PROPERTY(
    TARGET ${DTK_SEARCH_LIBNAME}
    APPEND PROPERTY LINKER_LANGUAGE CXX
    )
ENDIF()

ADD_SUBDIRECTORY(Discretization)
ADD_SUBDIRECTORY(Meshfree)
ADD_SUBDIRECTORY(Benchmarks)
