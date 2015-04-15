
#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Client)
APPEND_SET(HEADERS
  ${DIR}/DTK_Entity.hpp
  ${DIR}/DTK_EntityExtraData.hpp
  ${DIR}/DTK_EntityImpl.hpp
  ${DIR}/DTK_EntityIntegrationRule.hpp
  ${DIR}/DTK_EntityIterator.hpp
  ${DIR}/DTK_EntityLocalMap.hpp
  ${DIR}/DTK_Field.hpp
  ${DIR}/DTK_EntitySet.hpp
  ${DIR}/DTK_EntityShapeFunction.hpp
  ${DIR}/DTK_Types.hpp
  ) 

APPEND_SET(SOURCES
  ${DIR}/DTK_Entity.cpp
  ${DIR}/DTK_EntityIterator.cpp
  ${DIR}/DTK_EntityLocalMap.cpp
  ${DIR}/DTK_EntitySet.cpp
  ${DIR}/DTK_EntityShapeFunction.cpp
  )

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Operator)
APPEND_SET(HEADERS
  ${DIR}/DTK_BasicEntityPredicates.hpp
  ${DIR}/DTK_FieldMultiVector.hpp
  ${DIR}/DTK_FieldMultiVector_impl.hpp
  ${DIR}/DTK_FunctionSpace.hpp
  ${DIR}/DTK_MapOperator.hpp
  ${DIR}/DTK_MapOperator_impl.hpp
  ) 

APPEND_SET(SOURCES
  ${DIR}/DTK_BasicEntityPredicates.cpp
  ${DIR}/DTK_FunctionSpace.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_interface
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
