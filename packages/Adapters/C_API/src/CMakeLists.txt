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
  DTK_C_API.h
  ) 

APPEND_SET(SOURCES
  DTK_POD_PointCloudEntity.cpp
  DTK_POD_PointCloudEntityImpl.cpp
  DTK_POD_PointCloudEntityIterator.cpp
  DTK_POD_PointCloudEntitySet.cpp
  DTK_POD_PointCloudLocalMap.cpp  
  DTK_C_API.cpp
  ) 

TRIBITS_ADD_LIBRARY(
  dtk_c_api
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
