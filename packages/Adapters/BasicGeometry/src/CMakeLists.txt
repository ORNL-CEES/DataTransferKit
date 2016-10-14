
#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  DTK_BasicEntitySet.hpp
  DTK_BasicGeometryEntity.hpp
  DTK_BasicGeometryManager.hpp
  DTK_EntityCenteredField.hpp
  DTK_EntityCenteredShapeFunction.hpp
  DTK_BasicGeometryEntityImpl.hpp
  DTK_BasicGeometryExtraData.hpp
  DTK_BasicGeometryLocalMap.hpp
  DTK_BoxGeometry.hpp
  DTK_BoxGeometryImpl.hpp
  DTK_CylinderGeometry.hpp
  DTK_CylinderGeometryImpl.hpp
  DTK_Point.hpp
  DTK_PointImpl.hpp
  )

APPEND_SET(SOURCES
  DTK_BasicEntitySet.cpp
  DTK_BasicGeometryEntity.cpp
  DTK_BasicGeometryManager.cpp
  DTK_EntityCenteredField.cpp
  DTK_EntityCenteredShapeFunction.cpp
  DTK_BasicGeometryExtraData.cpp
  DTK_BasicGeometryLocalMap.cpp
  DTK_BoxGeometry.cpp
  DTK_BoxGeometryImpl.cpp
  DTK_CylinderGeometry.cpp
  DTK_CylinderGeometryImpl.cpp
  DTK_Point.cpp
  DTK_PointImpl.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_basicgeometryadapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
