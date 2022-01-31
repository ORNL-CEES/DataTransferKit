INCLUDE(TribitsCreateClientTemplateHeaders)

#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRIBITS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})


# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_SET(HEADERS ${PACKAGE_BINARY_DIR}/${PACKAGE_NAME}_config.hpp)
APPEND_SET(HEADERS ${PACKAGE_BINARY_DIR}/${PACKAGE_NAME}_ETIHelperMacros.h)

# Explicitly instantiate classes.
IF (${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)

  # Generate ETI .cpp files for DataTransferKit::NearestNeighborOperator
  DTK_PROCESS_ALL_N_TEMPLATES(NEARESTNEIGHBOROPERATOR_OUTPUT_FILES
          "DTK_ETI_NT.tmpl" "NearestNeighborOperator" "NEARESTNEIGHBOROPERATOR"
    "${${PACKAGE_NAME}_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${NEARESTNEIGHBOROPERATOR_OUTPUT_FILES})

  # Generate ETI .cpp files for DataTransferKit::MovingLeastSquaresOperator
  DTK_PROCESS_ALL_N_TEMPLATES(MOVINGLEASTSQUARESOPERATOR_OUTPUT_FILES
          "DTK_ETI_NT.tmpl" "MovingLeastSquaresOperator" "MOVING_LEAST_SQUARES_OPERATOR"
    "${${PACKAGE_NAME}_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${MOVINGLEASTSQUARESOPERATOR_OUTPUT_FILES})

  # Generate ETI .cpp files for DataTransferKit::SplineOperator
  DTK_PROCESS_ALL_N_TEMPLATES(SPLINEOPERATOR_OUTPUT_FILES
          "DTK_ETI_NT.tmpl" "SplineOperator" "SPLINE_OPERATOR"
    "${${PACKAGE_NAME}_ETI_NODES}" TRUE)
  LIST(APPEND SOURCES ${SPLINEOPERATOR_OUTPUT_FILES})

ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_meshfree
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS dtk_utils dtk_search
  ADDED_LIB_TARGET_NAME_OUT DTK_MESHFREE_LIBNAME
  )

# We need to set the linker language explicitly here for CUDA builds.
SET_PROPERTY(
  TARGET ${DTK_MESHFREE_LIBNAME}
  APPEND PROPERTY LINKER_LANGUAGE CXX
  )
