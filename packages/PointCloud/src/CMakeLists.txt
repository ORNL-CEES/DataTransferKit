INCLUDE(TribitsLibraryMacros)

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
  DTK_BuhmannBasis.hpp
  DTK_BuhmannBasis_impl.hpp
  DTK_CenterDistributor.hpp
  DTK_CenterDistributor_impl.hpp
  DTK_CloudDomain.hpp
  DTK_CloudDomain_impl.hpp
  DTK_EuclideanDistance.hpp
  DTK_EuclideanDistance_impl.hpp
  DTK_LocalMLSProblem.hpp
  DTK_LocalMLSProblem_impl.hpp
  DTK_MeshFreeInterpolator.hpp
  DTK_MeshFreeInterpolatorFactory.hpp
  DTK_MeshFreeInterpolatorFactory_impl.hpp
  DTK_MovingLeastSquare.hpp
  DTK_MovingLeastSquare_impl.hpp
  DTK_PointCloudDummy.hpp
  DTK_PolynomialMatrix.hpp
  DTK_PolynomialMatrix_impl.hpp
  DTK_RadialBasisPolicy.hpp
  DTK_SplineInterpolationPairing.hpp
  DTK_SplineInterpolationPairing_impl.hpp
  DTK_SplineInterpolator.hpp
  DTK_SplineInterpolator_impl.hpp
  DTK_SplineOperatorA.hpp
  DTK_SplineOperatorA_impl.hpp
  DTK_SplineOperatorC.hpp
  DTK_SplineOperatorC_impl.hpp
  DTK_WendlandBasis.hpp
  DTK_WendlandBasis_impl.hpp
  DTK_WuBasis.hpp
  DTK_WuBasis_impl.hpp
  ) 

APPEND_SET(SOURCES
  DTK_PointCloudDummy.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  dtk_pointcloud
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/
)
