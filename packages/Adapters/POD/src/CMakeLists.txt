SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  DTK_POD_PointCloudEntity.hpp
  DTK_POD_PointCloudEntityImpl.hpp
  DTK_POD_PointCloudEntityIterator.hpp
  DTK_POD_PointCloudEntitySet.hpp
  DTK_POD_PointCloudLocalMap.hpp
  DTK_POD_Types.hpp
  ) 

APPEND_SET(SOURCES
  DTK_POD_PointCloudEntity.cpp
  DTK_POD_PointCloudEntityImpl.cpp
  DTK_POD_PointCloudEntityIterator.cpp
  DTK_POD_PointCloudEntitySet.cpp
  DTK_POD_PointCloudLocalMap.cpp  
  ) 

TRIBITS_ADD_LIBRARY(
  dtk_pod
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

TRIBITS_ADD_LIBRARY(
  dtk_c_api
  HEADERS DTK_C_API.h
  SOURCES DTK_C_API.cpp
  )
TARGET_LINK_LIBRARIES(dtk_c_api dtk_pod)

TRIBITS_ADD_LIBRARY(
  dtk_fortran_api
  HEADERS ""
  SOURCES DTK_Fortran_API.F90
  )
TARGET_LINK_LIBRARIES(dtk_fortran_api dtk_c_api)

