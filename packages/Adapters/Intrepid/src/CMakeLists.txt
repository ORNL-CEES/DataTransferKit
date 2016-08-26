
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
  DTK_IntrepidBasisFactory.hpp
  DTK_IntrepidCell.hpp
  DTK_IntrepidCellLocalMap.hpp
  DTK_IntrepidIntegrationRule.hpp
  DTK_IntrepidShapeFunction.hpp
  DTK_IntrepidSideCell.hpp
  DTK_NewtonSolver.hpp
  DTK_NewtonSolver_impl.hpp
  DTK_NonlinearProblemTraits.hpp
  DTK_ProjectionPrimitiveNonlinearProblems.hpp
  DTK_ProjectionPrimitives.hpp
  ) 

APPEND_SET(SOURCES
  DTK_IntrepidBasisFactory.cpp
  DTK_IntrepidCell.cpp
  DTK_IntrepidCellLocalMap.cpp
  DTK_IntrepidIntegrationRule.cpp
  DTK_IntrepidShapeFunction.cpp  
  DTK_IntrepidSideCell.cpp
  DTK_ProjectionPrimitiveNonlinearProblems.cpp
  DTK_ProjectionPrimitives.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_intrepidadapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
